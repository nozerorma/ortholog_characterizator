# Positional outputs consumed by PhyloPhere

VARIABILITY receives TRANSLATION.out.prot, published under PROT/bmge. Historical
staging names PROT_RAW and PROT_VAR_RAW do not change the coordinate space.

- HMMSCAN_BATCH scans the full-length human protein: the human CDS of the unfiltered codon alignment
  (`cds_dir`), translated with gaps removed (`full_length_query.py`), with the Pfam gathering thresholds
  (`hmmscan --cut_ga`). The BMGE-trimmed human sequence is not scanned, because trimming removes interior
  columns and the trimmed sequence joins stretches that are not contiguous in the protein.
- MAP_DOMAIN_VARIABILITY retains ali_start/ali_end as hmmscan query coordinates on the full-length protein
  and exports the trimmed MSA columns as alignment_positions, with coordinate_space=CAAS_0based. A hit
  residue is projected residue -> codon column of the unfiltered alignment (ori_codon_col) -> prot_ali_col
  through the MAP table; residues in removed columns are counted in n_removed and are not exported. Query
  length, endpoints, MAP width against the codon alignment, and entropy coverage against the retained MAP
  columns are checked, and any mismatch stops the export.
- The reporting criterion is recorded: MERGE_DOMTBLOUT keeps the `# Option settings` line of hmmscan, and the
  exporter requires `--cut_ga` in it, or an explicit `--i_evalue_threshold` (independent E-value) for a scan
  made without it.
- AGGREGATE_UCR retains its original 1-based MSA position and declares
  coordinate_space=trimmed_alignment_1based. Its selection joins use that same
  trimmed index; Pfam joins use the exported mask, not query endpoints.
- The query producer, the Pfam exporter and the UCR exporter include their Python source hash in the task
  script, so changing one invalidates its Nextflow cache. A preserved domtblout can be supplied with
  hmmscan_domtblout to re-export existing scan results without scanning again; it must come from the same
  full-length queries and `--cut_ga`. `-entry PFAM_ONLY` reruns only the Pfam scan and export of a finished run.

These changes are the producer/consumer contract, not a post-processing patch
on lists. Old tables without coordinate_space remain historical artifacts and
must be re-exported by the corrected producer before being treated as current.

Regression tests cover native exports with gaps/X, disagreement between original
and trimmed MAP coordinates, direct producer-to-PhyloPhere consumption in the
paired checkout, UCR/selection joins, and UCR/Pfam mask overlaps.
