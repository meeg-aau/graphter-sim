/*
 * Pangenome graph of the pangenome-group target genomes
 */

// Concatenate genomes with PanSN-spec sequence names: <genome>#1#<contig>
process PANSN_FASTA {
    label 'pangenome'
    conda "${projectDir}/envs/pggb.yml"

    input:
    path genomes, stageAs: 'genomes/*'
    path names

    output:
    path 'pangenome.fa.gz',     emit: fasta
    path 'pangenome.fa.gz.fai', emit: fai
    path 'pangenome.fa.gz.gzi', emit: gzi

    script:
    """
    n=\$(grep -c . ${names} || true)
    if [[ \$n -lt 2 ]]; then
        echo "The pangenome needs at least 2 genomes (got \$n); adjust the split" >&2
        exit 1
    fi
    while read -r g; do
        sed "s/^>\\([^[:space:]]*\\).*/>\${g}#1#\\1/" "genomes/\${g}.fna"
    done < ${names} | bgzip -@ ${task.cpus} -c > pangenome.fa.gz
    samtools faidx pangenome.fa.gz
    """

    stub:
    """
    touch pangenome.fa.gz pangenome.fa.gz.fai pangenome.fa.gz.gzi
    """
}

process PGGB {
    label 'pangenome'
    conda "${projectDir}/envs/pggb.yml"
    publishDir "${params.outdir}/pangenome", mode: params.publish_mode

    input:
    path fasta
    path fai
    path gzi
    path names

    output:
    path 'pangenome.gfa', emit: gfa
    path 'pggb_out/*.{og,html,png,log,yml}', optional: true

    script:
    """
    n=\$(grep -c . ${names})
    pggb -i ${fasta} -o pggb_out -t ${task.cpus} -n \$n ${params.pggb_args}
    cp \$(ls pggb_out/*.final.gfa | head -n 1) pangenome.gfa
    """

    stub:
    """
    touch pangenome.gfa
    """
}

process MINIGRAPH {
    label 'pangenome'
    conda "${projectDir}/envs/minigraph.yml"
    publishDir "${params.outdir}/pangenome", mode: params.publish_mode

    input:
    path genomes, stageAs: 'genomes/*'
    path names

    output:
    path 'pangenome.gfa', emit: gfa

    script:
    """
    # Graph is seeded with the first genome in the list, then others are added
    mkdir pansn
    while read -r g; do
        sed "s/^>\\([^[:space:]]*\\).*/>\${g}#1#\\1/" "genomes/\${g}.fna" > "pansn/\${g}.fna"
    done < ${names}
    minigraph ${params.minigraph_args} -t ${task.cpus} \\
        \$(sed 's/^\\(.*\\)\$/pansn\\/\\1.fna/' ${names}) > pangenome.gfa
    """

    stub:
    """
    touch pangenome.gfa
    """
}
