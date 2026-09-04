# ONTAnalysis project structure and development guide

This document is the source of truth for future work by Qwen, Codex, and human
contributors. If code, README text, or an output directory disagrees with this
guide, either bring it into compliance or update this guide explicitly as part
of the same reviewed change.

## 1. Project scope

This project has one responsibility: turn the existing `Pa2_gDNA-2` ONT
modified-base reads into validated, reusable ONT results.

Required final products are:

1. A coordinate-sorted, indexed GRCm38p6 ONT modBAM.
2. Alignment and sequencing-coverage QC.
3. Separate 5mC and 5hmC results, with CpG as the primary context.
4. Coverage-weighted global, chromosome, and genomic-window summaries.
5. Single-sample QC, genome-wide, regional, and genomic-feature visualizations,
   with exact plot-data tables and explicit interpretation boundaries.
6. Sufficient manifests, logs, counts, and intermediate results to reproduce
   and audit every reported value.

SixBase sample-data ingestion, DMR/DhMR selection, and ONT-versus-SixBase
comparisons are out of scope. Stage 05 may reuse the exact mm10 annotation
assets and compatible figure families from SixBase. A separate project may
consume the final BAM, BAI, coverage tables, and methylation tables produced
here.

## 2. Scientific invariants

- The original instrument output is immutable.
- The primary read input is the 74 pass BAM files, not FASTQ. BAM retains the
  `MM`, `ML`, and `MN` modified-base tags.
- The reference is standard mouse GRCm38p6 without synthetic-control contigs.
- The existing basecalls must not be repeated unless a separate, explicitly
  approved re-basecalling analysis is created.
- 5mC and 5hmC must remain separate in commands, counts, tables, and reports.
- Report underlying call counts and denominators alongside percentages.
- The primary global estimate is coverage-weighted:

  ```text
  modification fraction = sum(modified calls) / sum(valid cytosine calls)
  ```

- Sequencing coverage and valid modification-call coverage are distinct
  quantities and must never be treated as interchangeable.

## 3. Storage boundary

### Primary drive: code and tools

Root: `/home/kyle/Desktop/Adam/ONTAnalysis`

```text
ONTAnalysis/
├── README.md
├── bin/                     Stable executable links
├── config/
│   └── project.env          Canonical paths and resource settings
├── docs/
│   ├── PROJECT_STRUCTURE.md This source-of-truth document
│   ├── workflow.md          Short workflow overview
│   └── audit_*.md           Dated audit records
├── scripts/
│   ├── 01_align_modbam.sh
│   ├── 02_finalize_alignment.sh
│   ├── 03_ont_qc_and_coverage.sh
│   ├── 04_methylation_analysis.sh
│   ├── 05-1_methylation_exploration.sh
│   ├── 05-2_methylation_exploration.R
│   ├── create_data_layout.sh
│   ├── stage_tools.sh
│   └── lib/
│       ├── pipeline_common.sh
│       ├── summarize_mosdepth.py
│       └── summarize_cpg_bedmethyl.py
└── tools/                   Pinned, versioned installations and archives
```

Do not place BAM, BED, bedMethyl, coverage, reference, or other large generated
data in this tree.

### 20 TB drive: all data and outputs

Root: `/mnt/20TB_A/Adam/ONTAnalysis`

```text
/mnt/20TB_A/Adam/ONTAnalysis/
├── Wei_P2a_modSUP/          Immutable original instrument output
├── inputs/
│   ├── reads/
│   │   └── Pa2_gDNA-2 -> ../../Wei_P2a_modSUP/.../run_directory
│   └── references/
│       └── GRCm38p6/
│           ├── GRCm38p6.fa
│           ├── GRCm38p6.fa.fai
│           └── GRCm38p6.lr-hq.v2.31.mmi
└── outputs/
    └── Pa2_gDNA-2/
        ├── 00_manifests/
        ├── 01_alignment/
        │   └── aligned_unsorted_chunks/
        ├── 02_alignment_qc/
        ├── 03_ont_qc_coverage/
        ├── 04_methylation/
        ├── 05_methylation_exploration/
        ├── logs/
        └── tmp/
```

Older empty directories such as `03_modkit`, `04_methylation_summary`, and
`05_gap_analysis` are not part of the canonical structure. Existing data must
never be deleted merely to make the tree match this guide; cleanup requires a
separate review of the exact paths and contents.

## 4. Five-stage pipeline contract

Each stage has one clear input/output contract. A later stage must fail safely
when the preceding stage is incomplete.

### Stage 01 — restartable alignment

Script: `scripts/01_align_modbam.sh`

Inputs:

- 74 unaligned pass modBAM files.
- GRCm38p6 FASTA.
- Minimap2 `lr:hq` index.

Tools:

- Standalone Minimap2 builds the `.mmi` index.
- Dorado `aligner` uses its bundled Minimap2 implementation for alignment.
- Samtools performs FASTA indexing and BAM integrity checks.

Required behavior:

- Align one source BAM at a time with Dorado `--no-sort`.
- Preserve `MM`, `ML`, and `MN` tags.
- Write to a partial filename, fully decompress it, validate its modification
  tags, then rename it atomically.
- Record the validation method, alignment-record count, byte size, modification
  time, and validation time in a versioned sidecar completion marker.
- Skip an unchanged chunk after a restart only when its marker matches the BAM.
- Recoverably quarantine and regenerate any completed chunk that fails a full
  stream check; do not delete the damaged BAM automatically.
- Never use Dorado's whole-run internal coordinate sort for this dataset.

Output:

```text
01_alignment/aligned_unsorted_chunks/
└── <source-name>.aligned.unsorted.bam
```

Completion gate: exactly 74 valid aligned chunks.

### Stage 02 — final alignment BAM

Script: `scripts/02_finalize_alignment.sh`

Input: the exact 74 stage-01 aligned chunks.

Tool: Samtools.

Required behavior:

- Refuse missing, extra, truncated, or unexpected chunks.
- Fully decompress every chunk before concatenation and write the per-chunk
  byte size, alignment-record count, and status to a validation manifest.
- Verify `MM`, `ML`, and `MN` tags in every chunk.
- Concatenate and coordinate-sort once with bounded CPU and memory.
- Fully decompress the sorted partial BAM and require its alignment-record count
  to equal the sum of the 74 input counts before promotion.
- Write BAM and BAI through partial files and atomic renames.
- Produce a manifest and alignment QC summary.

Canonical outputs:

```text
02_alignment_qc/
├── Pa2_gDNA-2.aligned.sorted.bam
├── Pa2_gDNA-2.aligned.sorted.bam.bai
└── Pa2_gDNA-2.alignment_summary.txt
```

Completion gate: BAM passes `samtools quickcheck` and a complete `samtools view`
stream, its alignment-record count matches the input total, header is
`SO:coordinate`, BAI is readable, and modification tags remain present.

### Stage 03 — ONT alignment QC and sequencing coverage

Script: `scripts/03_ont_qc_and_coverage.sh`

Input: final stage-02 BAM and BAI.

Tools:

- Samtools 1.22.1 reports alignment counts and statistics.
- Mosdepth 0.3.14 performs one efficient depth pass. Its standalone release
  binary is pinned by SHA-256 in `scripts/stage_tools.sh`.
- BGZF and Tabix from the Samtools environment compress and index the derived
  genomic tables.
- `scripts/lib/summarize_mosdepth.py` uses only the Python standard library to
  validate intervals and calculate tables. Its tiny-BAM results were audited
  against `samtools depth -aa` before full-data use.

Required analyses:

- Mapping, primary, secondary, supplementary, and unmapped counts.
- Mapping-quality and alignment statistics.
- Genome-wide and per-chromosome depth.
- Coverage breadth at at least 1x, 5x, 10x, and 20x.
- Fixed-window coverage suitable for later regional comparisons.
- Zero-coverage and low-coverage interval summaries.
- Optional compressed coverage track for IGV, without retaining an unnecessary
  uncompressed per-base depth file.

Canonical outputs:

```text
03_ont_qc_coverage/
├── alignment_flagstat.txt
├── alignment_stats.txt
├── alignment_idxstats.tsv
├── coverage_summary.tsv
├── chromosome_coverage.tsv
├── window_coverage.tsv.gz
├── window_coverage.tsv.gz.tbi
├── per_base_coverage.bed.gz
├── per_base_coverage.bed.gz.csi
├── zero_low_coverage_intervals.bed.gz
├── zero_low_coverage_intervals.bed.gz.tbi
├── mosdepth.global.dist.txt
├── mosdepth.region.dist.txt
├── mosdepth.summary.txt
└── ont_qc_coverage_report.txt
```

`per_base_coverage.bed.gz` is not a literal one-row-per-base file. Mosdepth
stores contiguous bases having the same integer depth as zero-based,
half-open run-length intervals, including depth-zero intervals. This preserves
localized gaps and provides an indexed IGV/query track without retaining a
large uncompressed depth file. The 100 kb window table reports mean depth,
base counts and breadth at 1x, 5x, 10x, and 20x, and an explicit window status.

Coverage uses minimum MAPQ 0 and excludes SAM flag mask 1796 (unmapped,
secondary, QC-fail, and duplicate records). Supplementary alignments are
included under Mosdepth's documented default policy. These settings are kept
in `config/project.env` and recorded in the report and manifest.

Completion gate: outputs are nonempty, internally consistent with the BAM's
mapped-read count, and explicitly include zero-coverage bases where relevant.

### Stage 04 — complete methylation analysis

Script: `scripts/04_methylation_analysis.sh`

This is intentionally one stage rather than separate Modkit check, pileup, and
summary scripts. Those operations share the same BAM, reference, configuration,
and failure boundary.

Inputs:

- Final stage-02 BAM and BAI.
- GRCm38p6 FASTA and FAI.

Primary tool: Modkit 0.6.4. Samtools validates the BAM and BGZF/Tabix index the
genomic tables. `scripts/lib/summarize_cpg_bedmethyl.py` is a standard-library
Python helper within Stage 04; it validates every paired bedMethyl record,
counts reference CpG sites, and performs the weighted aggregations.

Required behavior:

1. Validate modified-base tags.
2. Record all-context and CpG-focused Modkit summaries.
3. Produce a CpG bedMethyl result retaining separate 5mC and 5hmC records.
4. Calculate coverage-weighted global levels from counts.
5. Calculate chromosome and fixed-window levels.
6. Report valid, canonical, filtered, failed, and no-call counts where Modkit
   exposes them.
7. Validate every reported percentage against its numerator and denominator.

Canonical outputs:

```text
04_methylation/
├── modkit_tag_check.txt
├── modkit_valid_mm_headers.tsv
├── modkit_modified_bases.tsv
├── modkit_all_context_summary.txt
├── modkit_cpg_summary.txt
├── cpg_5mc_5hmc.bed.gz
├── cpg_5mc_5hmc.bed.gz.tbi
├── global_methylation_levels.tsv
├── chromosome_methylation_levels.tsv
├── window_methylation_levels.tsv.gz
├── window_methylation_levels.tsv.gz.tbi
├── methylation_call_coverage.tsv
└── methylation_analysis_report.txt
```

The full CpG bedMethyl is retained. A full all-context bedMethyl is disabled by
default because it is large and is not required for the initial goal; the
all-context summary is retained.

The CpG pileup uses
`--modified-bases 5mC 5hmC --cpg --combine-strands`. Thus one CpG dyad has
exactly one `m` row and one `h` row,
and both rows share the same valid-call denominator. Primary alignments are
used. Because this BAM contains localized coverage above Modkit's 65,535-count
limit, `--high-depth --max-depth 60000` prevents counter saturation while
retaining up to 60,000 observations per genomic position. The lowest 10th
percentile of modification probabilities is filtered, matching Modkit's
standard default while making both choices explicit in `config/project.env`.

`modkit_all_context_summary.txt` comes from `modkit summary --matched-only` and
is informational. `modkit_cpg_summary.txt` comes from `modkit stats` applied to
the indexed CpG pileup. Final 5mC and 5hmC quantities are calculated only from
the validated CpG bedMethyl counts. This avoids treating unmapped-tag totals as
CpG observations.

For every CpG pair, Stage 04 verifies:

```text
valid calls = 5mC + 5hmC + canonical C
5mC percent = sum(5mC calls) / sum(valid calls)
5hmC percent = sum(5hmC calls) / sum(valid calls)
```

It also reports reference CpG sites with and without a valid call, plus failed
or filtered, no-call, deletion, and different-base observations. These
modification-call coverage values are distinct from Stage 03 sequencing depth.

Completion gate: separate CpG 5mC and 5hmC counts and weighted percentages are
present, denominators are nonzero, and spot calculations reproduce the report.

### Stage 05 — methylation exploration and visualization

Scripts: `scripts/05-1_methylation_exploration.sh` and
`scripts/05-2_methylation_exploration.R`

Inputs:

- Completed Stage 03 coverage tables and Stage 04 methylation tables.
- The indexed, paired 5mC/5hmC CpG bedMethyl.
- GRCm38p6 FAI.
- The GENCODE vM25, CpG-island, cCRE, and intergenic annotation assets used by
  the SixBase mm10 workflow; their exact paths and checksums are manifested.

Required behavior:

- Produce publication-quality PDF versions of every figure.
- Produce a navigable HTML QC dashboard that consolidates Stage 02 alignment,
  Stage 03 sequencing coverage, Stage 04 modification calls, Stage 05 figures,
  exact plot-data links, provenance, and the original MinKNOW run report.
- Write the exact data used by each figure as TSV or compressed TSV.
- Preserve coverage-weighted 5mC and 5hmC calculations and keep sequencing
  depth distinct from valid modification-call depth.
- Reproduce the single-sample-compatible SixBase global, coverage, feature,
  scatter, and regional figure families.
- Label genome-position plots as descriptive rather than statistical and never
  present extreme windows as DMRs.
- Explicitly record why correlation, PCA, differential methylation, DMR
  heatmaps, and enrichment are unavailable with one sample.
- Support `--check-only`, validated skip, recoverable `--force` archiving,
  restartable feature aggregation, atomic promotion, logs, and manifests.

Canonical outputs:

```text
05_methylation_exploration/
├── .stage05.complete
├── ont_methylation_qc_report.html
├── methylation_exploration_report.txt
├── R_session_info.txt
├── figures/                       Publication-quality PDF figures
└── tables/                        15 plot-data/provenance tables
```

The primary figure groups are global CpG composition and levels, CpG coverage
and site distributions, chromosome and 100 kb window summaries, sequencing
versus valid-call coverage, descriptive genome-wide tracks, extreme-window
profiles, genomic-feature methylation boxplots, CpG-island/shore/shelf coverage
and length plots, and a weighted feature-class summary.
`sixbase_figure_mapping.tsv` identifies directly reproduced,
single-sample analogue, and scientifically unavailable SixBase families.

Completion gate: all expected PDFs have valid file signatures, the HTML and
plain-text reports, all 15 tables, and the session record are nonempty, and the
mapping table records the six multi-sample or differential families as unavailable.

## 5. Script engineering rules

All pipeline stages must:

- Begin with `set -Eeuo pipefail`.
- Source `config/project.env`; do not duplicate machine paths in scripts.
- Use `scripts/lib/pipeline_common.sh` for timestamped logging, locking, exit
  status, and resource heartbeats.
- Log tool versions, input paths, output paths, parameters, start/end times, and
  success or failure.
- Write logs to `outputs/Pa2_gDNA-2/logs`, not the primary drive.
- Offer a non-destructive preflight or `--check-only` mode when practical.
- Use bounded threads and memory so the workstation remains responsive.
- Write long-running outputs to `*.partial.*`, validate them, and atomically
  rename them only after success.
- Refuse to silently overwrite a validated final result.
- Be restartable at the smallest practical expensive unit.
- Fail on missing or unexpected inputs instead of continuing with a partial
  sample.
- Avoid retaining duplicate uncompressed/compressed genomic files.

## 6. Naming and manifests

- Sample identifier: `Pa2_gDNA-2`.
- Reference identifier: `GRCm38p6`.
- Use `.bam` only for completed BAMs; interrupted files include `.partial`.
- Use BED conventions explicitly: zero-based, half-open coordinates.
- TSV files must have a header and document units in either column names or the
  associated report.
- Manifests belong in `00_manifests` and should record paths, file sizes, tool
  versions, commands/options, and relevant checksums.
- Do not encode unrecorded thresholds only in filenames. Store them in the log,
  manifest, and report.

## 7. Retention and cleanup

Keep:

- Immutable original run data.
- Final aligned BAM and BAI.
- Compressed/indexed CpG bedMethyl.
- Final QC, coverage, methylation tables, manifests, and logs.
- Stage-01 chunks until the final BAM, Modkit outputs, and backups have been
  validated; their later deletion requires explicit approval.

Temporary files may be removed only after their corresponding final file has
passed validation. Interrupted or legacy outputs are not automatically deleted.
First inventory their count, size, path, and recoverability, then request an
explicit cleanup decision.

## 8. Change and review workflow

For Qwen or any implementation agent:

1. Read this document and `config/project.env` before editing.
2. State the exact stage and output contract being changed.
3. Keep unrelated changes separate.
4. Run syntax checks and a small integration test; never use the full dataset as
   the first test of a new command.
5. Do not start a full long-running stage as part of a code-review request.
6. Provide the diff, test command, test output, and expected full-run command to
   the reviewer.

For Codex or another reviewer:

1. Check scientific assumptions, tool semantics, coordinate conventions, and
   count denominators—not only shell syntax.
2. Confirm all outputs remain on the 20 TB drive and code remains on the primary
   drive.
3. Confirm restart, logging, partial-file, and completion-gate behavior.
4. Reproduce a small integration test and verify modification tags survive.
5. Inspect the relevant log and output inventory after each completed stage
   before approving the next one.

## 9. Current implementation status

- Stage 01: completed for all 74 source BAMs. One corrupt aligned output was
  detected by a full-stream check on 2026-08-28, recoverably quarantined, and
  regenerated from its intact source. All 74 current chunks have versioned,
  full-stream validation markers.
- Stage 02: implemented and audited; input and final-output validation were
  strengthened after the CRC incident. Ready to rerun on the repaired chunks.
- Stage 03: implemented and audited on a tiny BAM against Samtools; full run
  awaits the completed Stage 02 BAM.
- Stage 04: implemented and audited end to end on a tiny joint 5mC/5hmC modBAM;
  the full run awaits the completed Stage 02 BAM.
- SixBase comparison: intentionally outside this project.
