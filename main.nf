#!/usr/bin/env nextflow
/*
 * graphter-sim
 *
 * Synthetic data for testing sequence-graph tools:
 *   1. Download genomes of a target species and a background community (MiGA)
 *   2. Split target genomes into: metagenome-only, pangenome-only, both
 *   3. Simulate metagenomes (target + background) and assemble them (metaSPAdes)
 *   4. Build a pangenome graph of the pangenome genomes (PGGB or minigraph)
 *   5. Call genomovars among all target genomes (MiGA derep_wf, ANI 99.5%)
 */

nextflow.enable.dsl = 2

include { MIGA_LIST_TARGET                         } from './modules/miga'
include { MIGA_GTDB_GET as GET_TARGET              } from './modules/miga'
include { MIGA_GTDB_GET as GET_BACKGROUND          } from './modules/miga'
include { MIGA_DEREP                               } from './modules/miga'
include { STAGE_GENOMES as STAGE_TARGET            } from './modules/local'
include { STAGE_GENOMES as STAGE_BACKGROUND        } from './modules/local'
include { DESIGN; GENOMOVAR_TABLE                  } from './modules/local'
include { SIMULATE_READS                           } from './modules/simulate'
include { METASPADES                               } from './modules/assemble'
include { PANSN_FASTA; PGGB; MINIGRAPH             } from './modules/pangenome'

def fastaGlob(String p) {
    def f = file(p)
    return f.isDirectory() ? "${f}/*.{fna,fa,fasta,fas,fna.gz,fa.gz,fasta.gz,fas.gz}" : p
}

def validateParams() {
    if (!params.species && !params.target_genomes)
        error "Provide the target species (--species 'Genus species') or --target_genomes"
    def tech = params.read_technology.toString().toLowerCase()
    if (!(tech in ['illumina', 'illumina+nanopore', 'illumina+pacbio']))
        error "--read_technology must be one of: illumina, illumina+nanopore, illumina+pacbio"
    if (!(params.pangenome_tool in ['pggb', 'minigraph']))
        error "--pangenome_tool must be 'pggb' or 'minigraph'"
    if ((params.n_metagenomes as int) < 1)
        error "--n_metagenomes must be >= 1"
    // `--x_args '--flag'` makes Nextflow read --flag as a separate param and
    // set x_args to true, which would then be passed to the tool as "true"
    if (params.spades_ec_timeout) {
        try { params.spades_ec_timeout as Duration }
        catch (Exception e) { error "--spades_ec_timeout must be a duration such as '4.h' or '90.min', or 0 for no limit" }
    }
    params.keySet().findAll { it.endsWith('_args') }.each { k ->
        if (params[k] instanceof Boolean)
            error "--${k} needs a value. If it starts with '-', write it as --${k}='--option ...'"
    }
}

workflow {
    validateParams()
    // Effective parameters (defaults + profiles + command line), next to the
    // trace/report files of the same run
    def info_dir = file("${params.outdir}/pipeline_info")
    info_dir.mkdirs()
    def run_info = [
        run_name     : workflow.runName,
        session_id   : workflow.sessionId.toString(),
        resumed      : workflow.resume,
        command_line : workflow.commandLine,
        profile      : workflow.profile,
        nextflow     : nextflow.version.toString(),
        pipeline     : [version: workflow.manifest.version, commit: workflow.commitId ?: '', revision: workflow.revision ?: ''],
        params       : params.sort()
    ]
    file("${info_dir}/params_${params.trace_timestamp}.json").text =
        groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(run_info)) + '\n'
    // Printed up front: an onComplete handler has no syntax that compiles under
    // both the legacy parser (<=25.10) and the strict one (default in 26).
    log.info """
    graphter-sim outputs
      Metagenome assembly graphs : ${params.outdir}/metagenomes/assemblies/*.gfa
      Pangenome graph            : ${params.outdir}/pangenome/pangenome.gfa
      Genomovars                 : ${params.outdir}/genomovars/genomovars.tsv
      Ground truth / design      : ${params.outdir}/design/
    """.stripIndent()
    no_file = file("${projectDir}/assets/NO_FILE")

    // ---------------------------------------------------------------------
    // Target genomes
    // ---------------------------------------------------------------------
    // Full species list, used to keep the target species out of the background
    def download_bg = !params.background_genomes && (params.n_background as int) > 0
    if (params.species && download_bg) {
        MIGA_LIST_TARGET(params.species)
        target_all = MIGA_LIST_TARGET.out.list
    } else {
        target_all = Channel.value(no_file)
    }

    if (params.target_genomes) {
        STAGE_TARGET(
            Channel.fromPath(fastaGlob(params.target_genomes), checkIfExists: true)
                .collect().map { ['target', it] })
        target_genomes = STAGE_TARGET.out.genomes.map { it[1] }
        target_names   = STAGE_TARGET.out.names
    } else {
        def sp = params.species.toString().trim()
        GET_TARGET(
            Channel.of(['target', sp.startsWith('s__') ? sp : "s__${sp}",
                        params.max_target_genomes, false]),
            no_file)
        target_genomes = GET_TARGET.out.genomes.map { it[1] }
        target_names   = GET_TARGET.out.names
    }

    // ---------------------------------------------------------------------
    // Background community (one reference genome per random species)
    // ---------------------------------------------------------------------
    if (params.background_genomes) {
        STAGE_BACKGROUND(
            Channel.fromPath(fastaGlob(params.background_genomes), checkIfExists: true)
                .collect().map { ['background', it] })
        bg_genomes = STAGE_BACKGROUND.out.genomes.map { it[1] }
        bg_names   = STAGE_BACKGROUND.out.names
    } else if (download_bg) {
        if (!params.species)
            log.warn "No --species given: background genomes may include the target species"
        GET_BACKGROUND(
            Channel.of(['background', params.background_taxon, params.n_background, true]),
            target_all)
        bg_genomes = GET_BACKGROUND.out.genomes.map { it[1] }
        bg_names   = GET_BACKGROUND.out.names
    } else {
        bg_genomes = Channel.value([])
        bg_names   = Channel.value(no_file)
    }

    // ---------------------------------------------------------------------
    // Split target genomes and design the metagenomes
    // ---------------------------------------------------------------------
    DESIGN(target_names, bg_names)

    all_genomes = target_genomes.mix(bg_genomes).flatten().collect()

    // ---------------------------------------------------------------------
    // Metagenomes: simulate reads and assemble
    // ---------------------------------------------------------------------
    samples = DESIGN.out.design.flatMap { f ->
        (1..(params.n_metagenomes as int)).collect { k -> [String.format('mg%02d', k), f] }
    }
    SIMULATE_READS(samples, all_genomes)
    METASPADES(SIMULATE_READS.out.reads)

    // ---------------------------------------------------------------------
    // Pangenome graph
    // ---------------------------------------------------------------------
    if (params.pangenome_tool == 'pggb') {
        PANSN_FASTA(target_genomes, DESIGN.out.pangenome)
        PGGB(PANSN_FASTA.out.fasta, PANSN_FASTA.out.fai, PANSN_FASTA.out.gzi,
             DESIGN.out.pangenome)
    } else {
        MINIGRAPH(target_genomes, DESIGN.out.pangenome)
    }

    // ---------------------------------------------------------------------
    // Genomovars (all target genomes)
    // ---------------------------------------------------------------------
    MIGA_DEREP(target_genomes)
    GENOMOVAR_TABLE(MIGA_DEREP.out.clades, DESIGN.out.groups)

}
