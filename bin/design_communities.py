#!/usr/bin/env python3
"""
Partition the target genomes into three groups and design the composition of
the synthetic metagenomes.

Groups of target genomes:
  mg_only  genomes used only to simulate metagenomes
  pg_only  genomes used only to build the pangenome
  both     genomes used in both
  unused   leftover genomes (only when the split does not cover all genomes)

Outputs (in the current directory):
  target_groups.tsv      genome <TAB> group
  metagenome_design.tsv  sample, genome, source, group, depth
  pangenome_genomes.txt  genomes for the pangenome (pg_only + both)
"""

import argparse
import math
import random
import sys


def read_names(path):
    if path is None:
        return []
    with open(path) as fh:
        names = [ln.strip() for ln in fh if ln.strip() and not ln.startswith('#')]
    dups = {n for n in names if names.count(n) > 1}
    if dups:
        sys.exit(f'Duplicated genome names in {path}: {", ".join(sorted(dups))}')
    return sorted(names)


def split_counts(n, weights):
    """
    Number of genomes per group. If the weights add up to <= 1 they are
    fractions of n (largest-remainder rounding); otherwise they are counts.
    """
    total_w = sum(weights)
    if any(w < 0 for w in weights):
        sys.exit('Split values must be non-negative')
    if total_w <= 1 + 1e-6:
        total = n if abs(total_w - 1) < 1e-6 else int(round(total_w * n))
        exact = [w / total_w * total if total_w > 0 else 0 for w in weights]
        counts = [int(math.floor(x)) for x in exact]
        rema = sorted(range(len(weights)), key=lambda i: exact[i] - counts[i],
                      reverse=True)
        for i in rema[:total - sum(counts)]:
            counts[i] += 1
        return counts
    counts = [int(round(w)) for w in weights]
    if sum(counts) > n:
        sys.exit(f'Requested {sum(counts)} target genomes in the split, '
                 f'but only {n} are available')
    return counts


def lognormal_with_mean(rng, mean, sdlog):
    mu = math.log(mean) - sdlog ** 2 / 2
    return rng.lognormvariate(mu, sdlog)


def assign_targets(rng, pool, n_samples, per_sample):
    """
    Every genome in the pool is placed in at least one sample (round-robin on
    a shuffled pool), then samples are filled up to per_sample genomes.
    """
    pool = pool[:]
    rng.shuffle(pool)
    samples = [[] for _ in range(n_samples)]
    for i, g in enumerate(pool):
        samples[i % n_samples].append(g)
    if per_sample <= 0:
        per_sample = math.ceil(len(pool) / n_samples) if pool else 0
    per_sample = min(per_sample, len(pool))
    for s in samples:
        others = [g for g in pool if g not in s]
        rng.shuffle(others)
        while len(s) < per_sample and others:
            s.append(others.pop())
    return samples


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--target', required=True, help='File with target genome names')
    ap.add_argument('--background', help='File with background genome names')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--split-mg-only', type=float, default=1 / 3)
    ap.add_argument('--split-pg-only', type=float, default=1 / 3)
    ap.add_argument('--split-both', type=float, default=1 / 3)
    ap.add_argument('--n-metagenomes', type=int, default=3)
    ap.add_argument('--target-depth', type=float, default=50)
    ap.add_argument('--target-per-metagenome', type=int, default=0)
    ap.add_argument('--target-abundance', choices=['random', 'even'], default='random')
    ap.add_argument('--target-abundance-sdlog', type=float, default=1.0)
    ap.add_argument('--background-per-metagenome', type=int, default=0)
    ap.add_argument('--background-depth-mean', type=float, default=10)
    ap.add_argument('--background-depth-sdlog', type=float, default=1.0)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    target = read_names(args.target)
    background = read_names(args.background)
    clash = set(target) & set(background)
    if clash:
        sys.exit('Genomes both in target and background: ' + ', '.join(sorted(clash)))
    if not target:
        sys.exit('No target genomes available')
    if args.n_metagenomes < 1:
        sys.exit('--n-metagenomes must be >= 1')

    # Partition target genomes
    shuffled = target[:]
    rng.shuffle(shuffled)
    n_mg, n_pg, n_both = split_counts(
        len(target), [args.split_mg_only, args.split_pg_only, args.split_both])
    groups = {}
    i = 0
    for name, n in (('mg_only', n_mg), ('pg_only', n_pg), ('both', n_both)):
        for g in shuffled[i:i + n]:
            groups[g] = name
        i += n
    for g in shuffled[i:]:
        groups[g] = 'unused'

    with open('target_groups.tsv', 'w') as fh:
        fh.write('genome\tgroup\n')
        for g in target:
            fh.write(f'{g}\t{groups[g]}\n')

    pg_pool = sorted(g for g in target if groups[g] in ('pg_only', 'both'))
    with open('pangenome_genomes.txt', 'w') as fh:
        for g in pg_pool:
            fh.write(g + '\n')

    # Metagenome composition
    mg_pool = sorted(g for g in target if groups[g] in ('mg_only', 'both'))
    if not mg_pool:
        print('Warning: no target genomes assigned to metagenomes', file=sys.stderr)
    t_samples = assign_targets(rng, mg_pool, args.n_metagenomes,
                               args.target_per_metagenome)
    n_bg = args.background_per_metagenome
    if n_bg <= 0:
        n_bg = math.ceil(len(background) / 2)
    n_bg = min(n_bg, len(background))

    with open('metagenome_design.tsv', 'w') as fh:
        fh.write('sample\tgenome\tsource\tgroup\tdepth\n')
        for k in range(args.n_metagenomes):
            sample = f'mg{k + 1:02d}'
            tg = sorted(t_samples[k])
            if args.target_abundance == 'even':
                w = [1.0] * len(tg)
            else:
                w = [rng.lognormvariate(0, args.target_abundance_sdlog) for _ in tg]
            for g, wi in zip(tg, w):
                depth = args.target_depth * wi / sum(w)
                fh.write(f'{sample}\t{g}\ttarget\t{groups[g]}\t{depth:.4f}\n')
            for g in sorted(rng.sample(background, n_bg)):
                depth = lognormal_with_mean(rng, args.background_depth_mean,
                                            args.background_depth_sdlog)
                fh.write(f'{sample}\t{g}\tbackground\tNA\t{depth:.4f}\n')


if __name__ == '__main__':
    main()
