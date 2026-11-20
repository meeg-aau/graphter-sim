#!/usr/bin/env python3
"""
Convert MiGA derep_wf's genomospecies.tsv (clades at the requested ANI
threshold) into a long table of genomes per genomovar, annotated with the
target-genome group (mg_only, pg_only, both, unused).

Usage: genomovars_table.py genomospecies.tsv target_groups.tsv > genomovars.tsv
"""

import sys


def main(clades_f, groups_f):
    groups = {}
    with open(groups_f) as fh:
        next(fh)
        for ln in fh:
            g, grp = ln.rstrip('\n').split('\t')
            groups[g] = grp

    print('genomovar\trepresentative\tgenome\tis_representative\tgroup')
    seen = set()
    with open(clades_f) as fh:
        next(fh)
        for k, ln in enumerate(fh, start=1):
            f = ln.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            gvar = f'gvar_{k}'
            rep = f[1]
            for g in f[2].split(','):
                seen.add(g)
                print('\t'.join([gvar, rep, g, str(g == rep).lower(),
                                 groups.get(g, 'NA')]))
    # Genomes that MiGA did not place in any clade (e.g., failed or inactive)
    for g in sorted(set(groups) - seen):
        print('\t'.join(['NA', 'NA', g, 'false', groups[g]]))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
