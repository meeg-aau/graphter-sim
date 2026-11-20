/*
 * Metagenome assembly with metaSPAdes
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

    script:
    def tech  = params.read_technology.toString().toLowerCase()
    def lflag = tech.endsWith('nanopore') ? '--nanopore' : tech.endsWith('pacbio') ? '--pacbio' : ''
    def mem   = Math.max(1, task.memory.toGiga() as int)
    def ec_s  = params.spades_ec_timeout ? (params.spades_ec_timeout as Duration).toSeconds() : 0
    """
    lr=""
    if [[ -n "${lflag}" && -s "${long_reads}" ]]; then lr="${lflag} ${long_reads}"; fi
    metaspades.py -1 ${r1} -2 ${r2} \$lr -o spades \\
        -t ${task.cpus} -m ${mem} ${params.spades_args} &
    pid=\$!

    # Watchdog for read error correction (BayesHammer). On some read sets it
    # slows to near zero throughput and would run for days. If assembly has
    # not started after --spades_ec_timeout, stop metaSPAdes and fail with a
    # hint. Exit 124 is not in the retry list, so the task is not repeated.
    kill_tree() { local c; for c in \$(pgrep -P "\$1"); do kill_tree "\$c"; done; kill -TERM "\$1" 2>/dev/null || true; }
    ec_limit=${ec_s}
    if (( ec_limit > 0 )); then
        while kill -0 \$pid 2>/dev/null; do
            if grep -q '^===== Assembling started' spades/spades.log 2>/dev/null; then break; fi
            if (( SECONDS >= ec_limit )); then
                kill_tree \$pid
                wait \$pid 2>/dev/null || true
                cp spades/spades.log ${sample}.spades.log 2>/dev/null || true
                echo "ERROR: metaSPAdes read error correction did not finish within ${params.spades_ec_timeout} for ${sample}." >&2
                echo "BayesHammer can stall on some read sets (see README). Either raise --spades_ec_timeout" >&2
                echo "or skip error correction with: -resume --spades_args=--only-assembler" >&2
                exit 124
            fi
            sleep 30
        done
    fi
    wait \$pid

    cp spades/scaffolds.fasta ${sample}.scaffolds.fasta
    cp spades/contigs.fasta   ${sample}.contigs.fasta
    cp spades/assembly_graph_with_scaffolds.gfa ${sample}.assembly_graph_with_scaffolds.gfa
    cp spades/spades.log ${sample}.spades.log

    # Published files are copied above; the rest of spades/ is only needed for debugging
    if [[ "${params.keep_spades_dir}" != "true" ]]; then rm -rf spades; fi
    """

    stub:
    """
    touch ${sample}.scaffolds.fasta ${sample}.contigs.fasta ${sample}.spades.log
    touch ${sample}.assembly_graph_with_scaffolds.gfa
    """
}
