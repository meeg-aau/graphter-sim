#!/usr/bin/env python3
"""
Resolve IUPAC ambiguity codes in the sequence lines of a FASTQ stream.

Each code (R Y K M S W B D H V, either case) becomes one of the bases it
allows, chosen at random; case is kept. ACGTN are left alone and any other
character becomes N. Headers, '+' lines and qualities are copied as they are,
so read lengths do not change. Expects 4-line records (one sequence line).

The random choices are seeded from --seed and --key, so the same input gives
the same output.

Usage:
  resolve_iupac.py --seed SEED --key 'SAMPLE GENOME R1' --counts FILE \\
      < in.fastq > out.fastq
"""

import argparse
import hashlib
import random
import sys

IUPAC = {'R': 'AG', 'Y': 'CT', 'K': 'GT', 'M': 'AC', 'S': 'CG', 'W': 'AT',
         'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG'}
CHOICES = {}
for code, bases in IUPAC.items():
    CHOICES[ord(code)] = bases.encode()
    CHOICES[ord(code.lower())] = bases.lower().encode()
VALID = b'ACGTNacgtn'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[1])
    ap.add_argument('--seed', required=True, help='Base seed (e.g. --seed of the run)')
    ap.add_argument('--key', required=True, help='Text mixed into the seed (sample, genome, file)')
    ap.add_argument('--counts', required=True, help='File to write the number of replaced characters to')
    args = ap.parse_args()

    digest = hashlib.sha256(f'{args.seed}\t{args.key}'.encode()).digest()
    rng = random.Random(int.from_bytes(digest[:8], 'big'))

    n = 0
    out = sys.stdout.buffer
    for i, line in enumerate(sys.stdin.buffer):
        # Fast path: translate() in C drops valid bases, so most lines skip the loop
        if i % 4 == 1 and line.rstrip(b'\r\n').translate(None, VALID):
            seq = line.rstrip(b'\r\n')
            s = bytearray(seq)
            for j, c in enumerate(s):
                if c not in VALID:
                    bases = CHOICES.get(c)
                    s[j] = rng.choice(bases) if bases else ord('N')
                    n += 1
            line = bytes(s) + line[len(seq):]
        out.write(line)
    out.flush()

    with open(args.counts, 'w') as fh:
        fh.write(f'{n}\n')


if __name__ == '__main__':
    main()
