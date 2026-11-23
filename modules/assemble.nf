/*
 * Metagenome assembly with metaSPAdes
 *
 * Reads are error-corrected with non-meta `spades.py --only-error-correction`
 * and then assembled with `metaspades.py --only-assembler`. In metagenomic
 * mode, BayesHammer pre-filters singleton k-mers with a counting quotient
 * filter that can overfill on deep metagenomes: correction then crawls, or
 * finishes after silently dropping genuine k-mers (see README). Non-meta
 * error correction does not use that filter.
 */
process METASPADES {
    tag "${sample}"
    label 'assembly'
    conda "${projectDir}/envs/spades.yml"
    publishDir "${params.outdir}/metagenomes/assemblies", mode: params.publish_mode

    input:
    tuple val(sample), path(r1), path(r2), path(long_reads)

    output:
    tuple val(sample), path("${sample}.scaffolds.fasta"), emit: scaffolds
    tuple val(sample), path("${sample}.assembly_graph_with_scaffolds.gfa"), emit: gfa
    path "${sample}.contigs.fasta"
    path "${sample}.spades.log"
    path "${sample}.spades_ec.log", optional: true

    script:
    def tech     = params.read_technology.toString().toLowerCase()
    def lflag    = tech.endsWith('nanopore') ? '--nanopore' : tech.endsWith('pacbio') ? '--pacbio' : ''
    def mem      = Math.max(1, task.memory.toGiga() as int)
    def ec_s     = params.spades_ec_timeout ? (params.spades_ec_timeout as Duration).toSeconds() : 0
    def run_ec   = !params.spades_args.toString().tokenize().contains('--only-assembler')
    def asm_only = run_ec ? '--only-assembler' : ''  // otherwise already in spades_args
    """
    kill_tree() { local c; for c in \$(pgrep -P "\$1"); do kill_tree "\$c"; done; kill -TERM "\$1" 2>/dev/null || true; }
    # First path listed under key \$1 (e.g. 'left reads') in SPAdes' corrected.yaml.
    # Keys and values may be quoted: - "left reads":\\n  - "/path/R1...cor.fastq.gz"
    yaml_reads() {
        awk -v key="\$1" '
            /^[ -]*"?[a-z ]+"?:/ { k = \$0; sub(/^[ -]*"?/, "", k); sub(/"?:.*/, "", k); inkey = (k == key); next }
            inkey && /^ *- / { p = \$0; sub(/^ *- */, "", p); gsub(/"/, "", p); print p; exit }
        ' ec/corrected/corrected.yaml
    }

    r1c=${r1}; r2c=${r2}; sr=""
    if ${run_ec}; then
        # Step 1: error correction, short reads only. Long reads and
        # --spades_args belong to the assembly step.
        spades.py --only-error-correction -1 ${r1} -2 ${r2} -o ec \\
            -t ${task.cpus} -m ${mem} &
        pid=\$!

        # Safety net: stop error correction after --spades_ec_timeout. kill_tree
        # also stops spades-hammer, which coreutils' timeout would orphan.
        # Exit 124 is not in the retry list, so the task is not repeated.
        ec_limit=${ec_s}
        if (( ec_limit > 0 )); then
            while kill -0 \$pid 2>/dev/null; do
                if (( SECONDS >= ec_limit )); then
                    kill_tree \$pid
                    wait \$pid 2>/dev/null || true
                    cp ec/spades.log ${sample}.spades_ec.log 2>/dev/null || true
                    echo "ERROR: read error correction did not finish within ${params.spades_ec_timeout} for ${sample}." >&2
                    echo "Either raise --spades_ec_timeout or skip error correction with:" >&2
                    echo "  -resume --spades_args=--only-assembler" >&2
                    exit 124
                fi
                sleep 5
            done
        fi
        wait \$pid
        cp ec/spades.log ${sample}.spades_ec.log

        r1c=\$(yaml_reads 'left reads')
        r2c=\$(yaml_reads 'right reads')
        su=\$(yaml_reads 'single reads')
        if [[ ! -s "\$r1c" || ! -s "\$r2c" ]]; then
            echo "ERROR: corrected paired reads not found in ec/corrected/corrected.yaml" >&2
            exit 1
        fi
        # Unpaired corrected reads (mates whose partner was discarded), if any
        if [[ -n "\$su" && -s "\$su" && -n "\$(gzip -cdf "\$su" | head -c 1)" ]]; then sr="-s \$su"; fi

        # BayesHammer occasionally writes invalid bytes (e.g. 0xff) into
        # corrected sequences. metaspades.py then either crashes on them
        # (UnicodeDecodeError) or passes them on silently. Replace every byte
        # other than ACGTN on sequence lines with N; headers and qualities are
        # untouched and read lengths do not change. LC_ALL=C matches bytes
        # >= 0x80 one by one. Output is uncompressed to save CPU.
        sanitise() (  # <in.fastq[.gz]> <out.fastq> <file for the count>
            set -o pipefail
            gzip -cdf "\$1" | LC_ALL=C awk -v cnt="\$3" '
                NR % 4 == 2 { n += gsub(/[^ACGTNacgtn]/, "N") }
                { print }
                END { print n + 0 > cnt }' > "\$2"
        )
        echo 0 > ec/unpaired.nbad
        sanitise "\$r1c" ec/R1.clean.fastq ec/R1.nbad & p1=\$!
        sanitise "\$r2c" ec/R2.clean.fastq ec/R2.nbad & p2=\$!
        pids="\$p1 \$p2"
        if [[ -n "\$sr" ]]; then
            sanitise "\$su" ec/unpaired.clean.fastq ec/unpaired.nbad & pids="\$pids \$!"
            sr="-s ec/unpaired.clean.fastq"
        fi
        for p in \$pids; do
            if ! wait \$p; then echo "ERROR: failed to sanitise the corrected reads" >&2; exit 1; fi
        done
        r1c=ec/R1.clean.fastq; r2c=ec/R2.clean.fastq
        a=\$(cat ec/R1.nbad); b=\$(cat ec/R2.nbad); c=\$(cat ec/unpaired.nbad)
        msg="graphter-sim: replaced \$((a + b + c)) invalid bases with N (R1=\$a, R2=\$b, unpaired=\$c)"
        echo "\$msg" >> ${sample}.spades_ec.log
        if (( a + b + c > 0 )); then echo "WARN: ${sample}: \$msg" >&2; fi
    fi

    # Step 2: assembly
    lr=""
    if [[ -n "${lflag}" && -s "${long_reads}" ]]; then lr="${lflag} ${long_reads}"; fi
    metaspades.py ${asm_only} -1 \$r1c -2 \$r2c \$sr \$lr -o spades \\
        -t ${task.cpus} -m ${mem} ${params.spades_args}

    cp spades/scaffolds.fasta ${sample}.scaffolds.fasta
    cp spades/contigs.fasta   ${sample}.contigs.fasta
    cp spades/assembly_graph_with_scaffolds.gfa ${sample}.assembly_graph_with_scaffolds.gfa
    cp spades/spades.log ${sample}.spades.log

    # Published files are copied above; ec/ and spades/ are only needed for debugging
    if [[ "${params.keep_spades_dir}" != "true" ]]; then rm -rf spades ec; fi
    """

    stub:
    """
    touch ${sample}.scaffolds.fasta ${sample}.contigs.fasta ${sample}.spades.log
    touch ${sample}.spades_ec.log ${sample}.assembly_graph_with_scaffolds.gfa
    """
}
