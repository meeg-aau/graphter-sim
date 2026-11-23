#!/usr/bin/env bash
#
# Simulate the reads of one metagenome, genome by genome, so that each genome
# reaches exactly the depth requested in the design table.
#
# Usage:
#   simulate_metagenome.sh SAMPLE DESIGN_TSV GENOME_DIR TECHNOLOGY \
#       ILLUMINA_MODEL READ_LENGTH SEED CPUS LONG_DEPTH_SCALE LONG_ARGS
#
# Outputs:
#   SAMPLE_R1.fastq.gz, SAMPLE_R2.fastq.gz       Illumina pairs (InSilicoSeq)
#   SAMPLE_long.fastq.gz                         Long reads (Badread; empty if not hybrid)
#   SAMPLE.abundance.tsv                         Realized composition
#
# IUPAC ambiguity codes in the genomes are copied into the reads by the
# simulators, but no sequencer outputs them, and BayesHammer turns them into
# invalid 0xff bytes. resolve_iupac.py replaces each with a random compatible
# base (seeded) in the reads only; the genomes are left unchanged. The number
# of replaced characters per genome is the iupac_resolved column.

set -euo pipefail

sample="$1"; design="$2"; gdir="$3"; tech="$4"; model="$5"
read_len="$6"; seed="$7"; cpus="$8"; long_scale="$9"; long_args="${10:-}"

case "$tech" in
  illumina) long="" ;;
  illumina+nanopore) long="nanopore" ;;
  illumina+pacbio) long="pacbio" ;;
  *) echo "Unsupported read technology: $tech" >&2; exit 1 ;;
esac

model=$(echo "$model" | tr '[:upper:]' '[:lower:]')
if [[ "$read_len" -le 0 ]]; then
  case "$model" in
    hiseq) read_len=125 ;;
    novaseq) read_len=150 ;;
    miseq*|nextseq) read_len=300 ;;
    *) echo "Unknown read length for model '$model'; set --illumina_read_length" >&2
       exit 1 ;;
  esac
fi

: > "${sample}_R1.fastq.gz"
: > "${sample}_R2.fastq.gz"
: > "${sample}_long.fastq.gz"   # stays empty unless hybrid
printf "sample\tgenome\tsource\tgroup\tdepth\tgenome_length\tread_length\tn_reads\tiupac_resolved\n" \
  > "${sample}.abundance.tsv"

k=0
tail -n +2 "$design" | awk -F'\t' -v s="$sample" '$1 == s' |
while IFS=$'\t' read -r smp genome source group depth; do
  k=$((k + 1))
  fasta="$gdir/${genome}.fna"
  [[ -s "$fasta" ]] || { echo "Missing genome: $fasta" >&2; exit 1; }
  glen=$(awk '!/^>/{ n += length($0) } END { print n + 0 }' "$fasta")
  # InSilicoSeq --n_reads is the total number of reads (both mates)
  n_reads=$(awk -v d="$depth" -v g="$glen" -v l="$read_len" \
    'BEGIN { n = int(d * g / l + 0.999999); if (n % 2) n++; if (n < 2) n = 2; print n }')
  echo "[$sample] $genome ($source): ${depth}X, $glen bp, $n_reads reads" >&2
  tmp="tmp_${k}"
  iss generate --draft "$fasta" --model "$model" --n_reads "$n_reads" \
    --cpus "$cpus" --seed "$((seed + k))" --output "$tmp" --quiet >&2
  n_iupac=0
  for r in 1 2; do
    f="${tmp}_R${r}.fastq.gz"; [[ -s "$f" ]] || f="${tmp}_R${r}.fastq"
    gzip -cdf "$f" |
      resolve_iupac.py --seed "$seed" --key "$sample $genome R$r" --counts "${tmp}_R${r}.iupac" |
      gzip -c >> "${sample}_R${r}.fastq.gz"
    n_iupac=$((n_iupac + $(cat "${tmp}_R${r}.iupac")))
  done
  rm -f "${tmp}"_*

  if [[ -n "$long" ]]; then
    ldepth=$(awk -v d="$depth" -v s="$long_scale" 'BEGIN { printf "%.4f", d * s }')
    # shellcheck disable=SC2086
    badread simulate --reference "$fasta" --quantity "${ldepth}x" --seed "$((seed + k))" $long_args \
      2> "badread_${k}.log" |
      resolve_iupac.py --seed "$seed" --key "$sample $genome long" --counts "long_${k}.iupac" |
      gzip -c >> "${sample}_long.fastq.gz"
    n_iupac=$((n_iupac + $(cat "long_${k}.iupac")))
    rm -f "long_${k}.iupac"
  fi

  if (( n_iupac > 0 )); then
    echo "[$sample] $genome: resolved $n_iupac IUPAC codes in the reads" >&2
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$smp" "$genome" "$source" "$group" \
    "$depth" "$glen" "$read_len" "$n_reads" "$n_iupac" >> "${sample}.abundance.tsv"
done
