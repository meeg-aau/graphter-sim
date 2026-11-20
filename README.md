# Graphter Simulator

A Nextflow workflow that generates synthetic genomic data for testing
sequence-graph tools. It produces metagenome assembly graphs and a pangenome
graph for a target species. It also produces a ground-truth table of genomovars,
which are clusters of genomes at 99.5% ANI.

```
target species (GTDB) ─┬─► split ─┬─ mg_only ─┐
                       │          ├─ both ────┼──► simulate reads ─► metaSPAdes ─► *.gfa, scaffolds
                       │          └─ pg_only ─┴──► pangenome graph ─► pangenome.gfa
                       │                      (both feeds both branches)
                       └─► MiGA derep_wf (ANI 99.5%) ─► genomovars.tsv
background species (GTDB) ───────► simulate reads
```

## Steps

1. **Target genomes.** `miga gtdb_get -T "s__<species>" --max-download N`
   downloads a random subset of up to `--max_target_genomes` genomes (default
   100) of the species.
2. **Background community.**
   `miga gtdb_get -T <background_taxon> --reference --max-download N`
   downloads up to `--n_background` species-representative genomes (default 50).
   Each one comes from a different random species in `--background_taxon`
   (default `d__Bacteria`). All genomes of the target species are passed to
   `--exclude`, so none of them end up in the background.
3. **Partition.** The target genomes are split at random into three groups:
   * `mg_only`: genomes used only in the metagenomes.
   * `pg_only`: genomes used only in the pangenome.
   * `both`: genomes used in both.

   The default split is ⅓ for each group.
4. **Metagenomes.** The workflow simulates `--n_metagenomes` metagenomes
   (default 3).
   * Each metagenome gets a random subset of the `mg_only` + `both` genomes.
     Every one of those genomes appears in at least one metagenome.
   * Each metagenome also gets a random subset of the background genomes.
   * The target species always adds up to `--target_depth` (default 50X) per
     metagenome. That depth is shared among the target genomes, at random
     (log-normal) or evenly.
   * Background genomes get log-normal depths.
   * Reads are simulated one genome at a time, so each genome gets exactly its
     intended depth.
5. **Assembly.** metaSPAdes assembles each metagenome. The outputs are
   `scaffolds.fasta` and `assembly_graph_with_scaffolds.gfa`. That graph is the
   simplified assembly graph plus the contig and scaffold paths (`P` and `J`
   lines); dropping those lines leaves exactly SPAdes'
   `assembly_graph_after_simplification.gfa`.
6. **Pangenome graph.** The `pg_only` + `both` genomes are named with PanSN
   (`genome#1#contig`) and built into a graph with PGGB (default) or minigraph.
7. **Genomovars.** All target genomes go through
   `miga derep_wf --threshold 99.5`. The resulting clades are written as a long
   table, annotated with each genome's group.

## Main outputs

| Path | Content |
|---|---|
| `metagenomes/assemblies/<mg>.assembly_graph_with_scaffolds.gfa` | Simplified metaSPAdes assembly graph (GFA), with scaffold paths |
| `metagenomes/assemblies/<mg>.scaffolds.fasta` | Scaffolds |
| `pangenome/pangenome.gfa` | Pangenome graph |
| `genomovars/genomovars.tsv` | Genomes per genomovar: `genomovar, representative, genome, is_representative, group` |
| `design/target_groups.tsv` | Group of each target genome (`mg_only`, `pg_only`, `both`, `unused`) |
| `design/metagenome_design.tsv` | Intended composition and depth of each metagenome |
| `metagenomes/reads/<mg>.abundance.tsv` | Realized composition: genome length, read length, number of reads |
| `metagenomes/reads/<mg>_R{1,2}.fastq.gz` | Simulated reads |
| `genomovars/derep/` | Raw MiGA `derep_wf` output (`genomospecies.tsv`, summaries) |
| `pipeline_info/params_<timestamp>.json` | Parameters actually used (defaults, profiles and command line combined), plus the run name, session ID, command line and Nextflow version. Its timestamp matches the run's `report_`, `trace_` and `timeline_` files. |

## Usage

```bash
nextflow run meeg-aau/graphter-sim -profile conda  # Nextflow >= 24.04 \
    --species "Escherichia coli" --outdir results
```

Useful options (defaults are in parentheses). Values for the `*_args` options
often start with `-`, so pass them with `=`, e.g.
`--spades_args='--only-assembler'`. Otherwise Nextflow reads the value as a
separate option.

| Option | Description |
|---|---|
| `--species` | Target species, either `Genus species` or `s__Genus species` |
| `--max_target_genomes` (100) | Maximum number of target genomes to download |
| `--n_background` (50), `--background_taxon` (`d__Bacteria`) | Background community |
| `--target_genomes`, `--background_genomes` | Directory or glob of local FastA files to use instead of downloading |
| `--split_mg_only`, `--split_pg_only`, `--split_both` (⅓ each) | Group sizes. Values that add up to ≤ 1 are fractions; anything larger is a number of genomes. Genomes left over go to the `unused` group. |
| `--n_metagenomes` (3) | Number of metagenomes |
| `--target_depth` (50) | Total depth of the target species per metagenome |
| `--target_per_metagenome` (0 = auto) | Target genomes per metagenome |
| `--target_abundance` (`random`) | How depth is shared among target genomes: `random` (log-normal, `--target_abundance_sdlog`) or `even` |
| `--background_per_metagenome` (0 = half), `--background_depth_mean` (10), `--background_depth_sdlog` (1) | Background composition |
| `--read_technology` (`illumina`) | `illumina`, `illumina+nanopore` or `illumina+pacbio`. metaSPAdes always needs short reads, so long reads are only used for hybrid assembly. |
| `--illumina_model` (`novaseq`) | InSilicoSeq error model: `hiseq`, `novaseq`, `miseq` or `nextseq` |
| `--long_read_depth_scale` (1.0), `--badread_*_args` | Long-read simulation (Badread) |
| `--pangenome_tool` (`pggb`), `--pggb_args` (`-p 95 -s 5000`) | Pangenome graph |
| `--genomovar_ani` (99.5) | ANI threshold for `derep_wf` |
| `--derep_min_qual` (`no`) | Keeps every genome. MiGA's own default (25) would drop low-quality genomes. |
| `--derep_args` (`--only-ani`) | Extra options for `derep_wf` |
| `--spades_ec_timeout` (`4.h`) | Maximum time for metaSPAdes' read error correction. If it runs longer, the task fails with a hint instead of running for days (see the caveat below). `0` means no limit. |
| `--keep_spades_dir` (false) | Keep metaSPAdes' working directory (`spades/`: corrected reads, per-k graphs, temporary files) in each assembly task's work directory, for debugging. By default it is deleted after the published files are copied, which saves about 1.5–3 GB per assembly. If error correction times out, it is always kept. |
| `--seed` (42) | Seed for the partition, the design and InSilicoSeq |
| `--max_cpus`, `--max_memory` (this host's CPUs and 95% of its RAM), `--max_time` (240 h) | Upper limit for any single task. Larger requests (e.g. metaSPAdes' 64 GB) are capped to these instead of failing. With `-profile slurm`, set them to the largest node. |

`-profile test` runs a small configuration with 12 target genomes, 6 background
genomes and 2 metagenomes at 10X.

### Software

* **MiGA** must be installed and on `PATH`. If `$MIGA_HOME` (or `--miga_home`)
  already contains an initialized MiGA, it is reused. Otherwise each MiGA task
  runs `miga init --auto` in its work directory. `derep_wf` needs MiGA's
  complete set of dependencies. You can also set `--miga_container`.
  * MiGA 1.4.3.1's bundled FastAAI needs **pyhmmer < 0.11** (so Python ≤ 3.12).
    From 0.11 on pyhmmer returns names as `str`, FastAAI's `.decode()` calls
    fail inside a bare `except`, and `derep_wf` exits 0 without clustering
    anything. `MIGA_DEREP` now fails in that case. A working conda setup:
    `python=3.12 pyhmmer=0.10.15`.
  * `MIGA_DEREP` runs with a copy of the MiGA configuration but no local
    databases, so `derep_wf` skips classifying each genome against Phyla_Lite.
    That step is irrelevant for ANI genomovars and can take hours on file
    systems with slow fsync. Your `MIGA_HOME` is not modified.
* **Other tools** (InSilicoSeq, Badread, SPAdes, PGGB, minigraph, samtools) come
  from the conda environments in `envs/`:
  * `-profile conda` (or `mamba`) creates the environments.
  * `-profile docker` (or `singularity`) builds containers from them with Wave.

### Notes and caveats

* MiGA itself picks the target subset at random (`--max-download`). That choice
  does not follow `--seed`, so the downloaded genomes change between runs. To
  get fully reproducible runs, keep the downloaded genomes and pass them back
  with `--target_genomes` and `--background_genomes`.
* `--only-ani` turns off the hAAI and AAI pre-screening in `derep_wf`. That is
  appropriate within a species (AAI ≥ 85%). Set `--derep_args ''` to restore
  MiGA's defaults.
* `n_reads` follows InSilicoSeq's convention: it is the total number of reads
  across both mates.

* metaSPAdes' read error correction (BayesHammer) can slow to near zero
  throughput on some simulated read sets. It keeps every thread busy while
  reading input at a few KB per minute. This was seen with SPAdes 4.3.0 and
  3.15.5 on `-profile test` (E. coli, seed 42, metagenome `mg01`). There, each
  half of the reads (one 72% GC background genome, and everything else) is
  corrected in minutes, but the combination stalls. `METASPADES` stops after
  `--spades_ec_timeout` (default 4 h) if assembly has not started, and fails
  with exit 124 without retrying. Then either raise the limit or rerun with
  `-resume --spades_args=--only-assembler` to skip error correction. The graphs
  will then keep more error-induced tips and bubbles.

## License

graphter-sim is distributed under the [Artistic License 2.0](LICENSE).
