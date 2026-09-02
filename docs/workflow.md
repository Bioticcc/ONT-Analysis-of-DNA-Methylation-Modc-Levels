# Workflow

The current dataset was already basecalled with a SUP model and the joint
all-context 5mC/5hmC modified-base model. The remaining stages are:

1. Validate the 74 pass modBAM chunks and generate a manifest.
2. Align each BAM independently to standard mouse GRCm38p6 with Dorado
   `--no-sort`, preserving `MM`, `ML`, and `MN` tags. Completed chunks are
   restartable and are never overwritten by a partial process.
3. Require all 74 chunks, concatenate them, then coordinate-sort and index the
   combined modBAM with controlled Samtools CPU and memory settings.
4. Run ONT alignment and sequencing-coverage QC with Samtools and one Mosdepth
   pass, retaining indexed exact run-length depth and 100 kb coverage windows.
5. Run the combined methylation stage: Modkit validation and summary, CpG
   pileup, and separate coverage-weighted 5mC/5hmC summaries. Reference CpG-site
   coverage and Modkit failed/no-call evidence are reported separately from raw
   sequencing depth.
6. Generate single-sample methylation exploration figures and their exact plot
   data, including the scientifically compatible figure families from the
   SixBase workflow and genomic-feature annotations used there.

POD5 and FASTQ files are retained as source material but are not primary input
to this run of the workflow.

Every stage writes a timestamped combined stdout/stderr log under the sample's
`logs` directory. Logs contain tool versions, command settings, resource
snapshots, per-chunk progress, and a final success/failure status.

See `PROJECT_STRUCTURE.md` for the authoritative stage contracts, canonical
outputs, development rules, and project boundary. SixBase comparisons are not
performed in this project; Stage 05 reuses its annotation sources and visual
design families without treating one ONT sample as a differential experiment.
