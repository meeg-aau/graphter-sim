/*
 * Read simulation: InSilicoSeq (Illumina) and Badread (long reads, hybrid mode)
 */
process SIMULATE_READS {
    tag "${sample}"
    label 'simulate'
    conda "${projectDir}/envs/simulate.yml"
    publishDir "${params.outdir}/metagenomes/reads", mode: params.publish_mode,
        saveAs: { fn -> (fn.endsWith('_long.fastq.gz') && params.read_technology == 'illumina') ? null : fn }

    input:
    tuple val(sample), path(design)
    path genomes, stageAs: 'genomes/*'

    output:
    tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz"),
          path("${sample}_long.fastq.gz"), emit: reads
    path "${sample}.abundance.tsv", emit: abundance

    script:
    def tech = params.read_technology.toString().toLowerCase()
    def long_args = tech.endsWith('nanopore') ? params.badread_nanopore_args :
                    tech.endsWith('pacbio')   ? params.badread_pacbio_args   : ''
    // Different seed for each metagenome
    def seed = (params.seed as long) * 1000 + (sample.replaceAll(/\D/, '') as long) * 100000
    """
    simulate_metagenome.sh '${sample}' '${design}' genomes '${tech}' \\
        '${params.illumina_model}' ${params.illumina_read_length} ${seed} \\
        ${task.cpus} ${params.long_read_depth_scale} '${long_args}'
    """

    stub:
    """
    touch ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz ${sample}_long.fastq.gz ${sample}.abundance.tsv
    """
}
