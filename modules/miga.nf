/*
 * Processes that use MiGA (https://github.com/bio-miga/miga)
 */

// Make sure MiGA is initialized. Uses --miga_home (or $MIGA_HOME / $HOME) if it
// already contains .miga_rc and .miga_daemon.json; otherwise it initializes a
// task-local MiGA home with default answers.
def migaSetup() {
    def home = params.miga_home ? "export MIGA_HOME='${params.miga_home}'" : ''
    def key  = params.ncbi_api_key ? "export NCBI_API_KEY='${params.ncbi_api_key}'" : ''
    """
    ${home}
    ${key}
    _mh="\${MIGA_HOME:-\$HOME}"
    if [[ ! -s "\$_mh/.miga_rc" || ! -s "\$_mh/.miga_daemon.json" ]]; then
        export MIGA_HOME="\$PWD/.miga_home"
        mkdir -p "\$MIGA_HOME"
        miga init --auto --no-optional --daemon-type bash < /dev/null \\
            > miga_init.log 2>&1 || { cat miga_init.log >&2; exit 1; }
    fi
    """.stripIndent()
}

def gtdbSpecies(String sp) {
    def s = sp.trim()
    return s.startsWith('s__') ? s : "s__${s}"
}

// Full list of genomes (dataset names) classified in the target species
process MIGA_LIST_TARGET {
    tag "${species}"
    label 'miga_download'

    input:
    val species

    output:
    path 'target_species_all.txt', emit: list

    script:
    def taxon = gtdbSpecies(species)
    """
    ${migaSetup()}
    miga new -P list_project -t genomes
    miga gtdb_get -P list_project -T '${taxon}' --dry -R target_species_all.txt
    [[ -s target_species_all.txt ]] || { echo "No genomes found in GTDB for '${taxon}'" >&2; exit 1; }
    """

    stub:
    """
    for i in \$(seq 1 ${params.max_target_genomes}); do echo "GCA_\${i}_1"; done > target_species_all.txt
    """
}

// Random subset of (up to) n genomes from a GTDB taxon
process MIGA_GTDB_GET {
    tag "${role}: ${taxon}"
    label 'miga_download'

    input:
    tuple val(role), val(taxon), val(n), val(reference_only)
    path exclude

    output:
    tuple val(role), path("genomes/*.fna"), emit: genomes
    path "${role}_genomes.txt", emit: names

    script:
    def ref  = reference_only ? '--reference' : ''
    def excl = exclude.name != 'NO_FILE' ? "--exclude '${exclude}'" : ''
    """
    ${migaSetup()}
    miga new -P project -t genomes
    miga gtdb_get -P project -T '${taxon}' --max-download ${n} ${ref} ${excl} \\
        ${params.miga_get_args}
    mkdir genomes
    for f in project/data/05.assembly/*.LargeContigs.fna*; do
        [[ -e "\$f" ]] || continue
        name=\$(basename "\$f" | sed -E 's/\\.LargeContigs\\.fna(\\.gz)?\$//')
        case "\$f" in
            *.gz) gzip -dc "\$f" > "genomes/\$name.fna" ;;
            *)    cp "\$f" "genomes/\$name.fna" ;;
        esac
    done
    ls genomes | sed 's/\\.fna\$//' > ${role}_genomes.txt
    [[ -s ${role}_genomes.txt ]] || { echo "No genomes downloaded for '${taxon}'" >&2; exit 1; }
    """

    stub:
    def nn = role == 'target' ? params.max_target_genomes : params.n_background
    def pf = role == 'target' ? 'GCA_' : 'BG_'
    """
    mkdir genomes
    for i in \$(seq 1 ${nn}); do
        printf ">c1\\nACGTACGTACGT\\n" > genomes/${pf}\${i}_1.fna
    done
    ls genomes | sed 's/\\.fna\$//' > ${role}_genomes.txt
    """
}

// Genomovars: dereplicate all target genomes at a high ANI threshold
process MIGA_DEREP {
    label 'miga_derep'
    publishDir "${params.outdir}/genomovars", mode: params.publish_mode

    input:
    path genomes, stageAs: 'genomes/*'

    output:
    path 'derep/genomospecies.tsv', emit: clades
    path 'derep/*.tsv'
    path 'derep/representatives.txt', optional: true

    script:
    """
    ${migaSetup()}
    # If Phyla_Lite is in MIGA_HOME, derep_wf's essential_genes step classifies
    # every genome with unknown domain against it (FastAAI hAAI). That is
    # irrelevant for ANI genomovars and very slow on file systems with slow
    # fsync. Run with the same MiGA configuration but no local databases. The
    # original MIGA_HOME is left untouched.
    _cfg="\${MIGA_HOME:-\$HOME}"
    export MIGA_HOME="\$PWD/.miga_home_nodb"
    mkdir -p "\$MIGA_HOME/.miga_db"
    cp "\$_cfg/.miga_rc" "\$_cfg/.miga_daemon.json" "\$MIGA_HOME/"
    printf '{"last_update":"%s","databases":{}}\\n' "\$(date -u '+%F %T +0000')" \\
        > "\$MIGA_HOME/.miga_db/_local_manif.json"

    miga derep_wf -o derep --threshold ${params.genomovar_ani} \\
        --min-qual ${params.derep_min_qual} \\
        -j ${task.cpus} -t 1 --threads-project ${task.cpus} \\
        ${params.derep_args} genomes/*.fna

    # derep_wf exits 0 even when every per-genome job failed, leaving a
    # genomospecies.tsv with only its header. Fail here instead of
    # publishing a genomovar table full of NA.
    n_in=\$(ls genomes/*.fna | wc -l)
    n_out=\$(awk -F'\\t' 'NR > 1 { n = split(\$3, m, ","); for (i = 1; i <= n; i++) seen[m[i]] = 1 }
                         END { print length(seen) }' derep/genomospecies.tsv)
    if [[ "\$n_out" -ne "\$n_in" ]]; then
        echo "ERROR: derep_wf clustered \$n_out of \$n_in genomes. Per-genome MiGA logs: \$PWD/derep/daemon/d/" >&2
        log=\$(grep -l -E 'Traceback|Error|failed' derep/daemon/d/*.log 2>/dev/null | head -1)
        [[ -n "\$log" ]] && { echo "--- tail of \$log ---" >&2; grep -v -E '^(DEBUG|OPT|CMD|\\s+from )' "\$log" | tail -15 >&2; }
        exit 1
    fi
    """

    stub:
    """
    mkdir derep
    printf "Clade\\tRepresentative\\tMembers\\n" > derep/genomospecies.tsv
    k=0
    for f in genomes/*.fna; do
        k=\$((k + 1)); n=\$(basename "\$f" .fna)
        printf "gsp_%s\\t%s\\t%s\\n" \$k \$n \$n >> derep/genomospecies.tsv
    done
    """
}
