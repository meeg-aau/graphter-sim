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
printf "sample\tgenome\tsource\tgroup\tdepth\tgenome_length\tread_length\tn_reads\n" \
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
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$smp" "$genome" "$source" "$group" \
    "$depth" "$glen" "$read_len" "$n_reads" >> "${sample}.abundance.tsv"

  echo "[$sample] $genome ($source): ${depth}X, $glen bp, $n_reads reads" >&2
  tmp="tmp_${k}"
  iss generate --draft "$fasta" --model "$model" --n_reads "$n_reads" \
    --cpus "$cpus" --seed "$((seed + k))" --output "$tmp" --quiet >&2
  for r in 1 2; do
    if [[ -s "${tmp}_R${r}.fastq.gz" ]]; then
      cat "${tmp}_R${r}.fastq.gz" >> "${sample}_R${r}.fastq.gz"
    else
      gzip -c "${tmp}_R${r}.fastq" >> "${sample}_R${r}.fastq.gz"
    fi
  done
  rm -f "${tmp}"_*

  if [[ -n "$long" ]]; then
    ldepth=$(awk -v d="$depth" -v s="$long_scale" 'BEGIN { printf "%.4f", d * s }')
    # shellcheck disable=SC2086
    badread simulate --reference "$fasta" --quantity "${ldepth}x" $long_args \
      2> "badread_${k}.log" | gzip -c >> "${sample}_long.fastq.gz"
  fi
done
