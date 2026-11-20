/*
 * Helper processes (no external bioinformatics tools)
 */

// Copy user-provided genomes as <sanitized_name>.fna (uncompressed)
process STAGE_GENOMES {
    tag "${role}"
    conda "${projectDir}/envs/python.yml"

    input:
    tuple val(role), path(files, stageAs: 'input/*')

    output:
    tuple val(role), path('genomes/*.fna'), emit: genomes
    path "${role}_genomes.txt", emit: names

    script:
    """
    mkdir genomes
    for f in input/*; do
        name=\$(basename "\$f" | sed -E 's/\\.gz\$//; s/\\.(fna|fasta|fa|fas|fsa|fnn)\$//' | sed 's/[^A-Za-z0-9_]/_/g')
        [[ -e "genomes/\$name.fna" ]] && { echo "Duplicated genome name: \$name" >&2; exit 1; }
        case "\$f" in
            *.gz) gzip -dc "\$f" > "genomes/\$name.fna" ;;
            *)    cp -L "\$f" "genomes/\$name.fna" ;;
        esac
    done
    ls genomes | sed 's/\\.fna\$//' > ${role}_genomes.txt
    """

    stub:
    """
    mkdir genomes
    for f in input/*; do touch "genomes/\$(basename "\$f" | sed -E 's/\\..*//').fna"; done
    ls genomes | sed 's/\\.fna\$//' > ${role}_genomes.txt
    """
}

process DESIGN {
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/design", mode: params.publish_mode

    input:
    path target_names
    path background_names

    output:
    path 'target_groups.tsv',     emit: groups
    path 'metagenome_design.tsv', emit: design
    path 'pangenome_genomes.txt', emit: pangenome

    script:
    """
    design_communities.py \\
        --target ${target_names} --background ${background_names} \\
        --seed ${params.seed} \\
        --split-mg-only ${params.split_mg_only} \\
        --split-pg-only ${params.split_pg_only} \\
        --split-both ${params.split_both} \\
        --n-metagenomes ${params.n_metagenomes} \\
        --target-depth ${params.target_depth} \\
        --target-per-metagenome ${params.target_per_metagenome} \\
        --target-abundance ${params.target_abundance} \\
        --target-abundance-sdlog ${params.target_abundance_sdlog} \\
        --background-per-metagenome ${params.background_per_metagenome} \\
        --background-depth-mean ${params.background_depth_mean} \\
        --background-depth-sdlog ${params.background_depth_sdlog}
    """
}

process GENOMOVAR_TABLE {
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/genomovars", mode: params.publish_mode

    input:
    path clades
    path groups

    output:
    path 'genomovars.tsv', emit: table

    script:
    """
    genomovars_table.py ${clades} ${groups} > genomovars.tsv
    """
}
