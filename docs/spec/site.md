# Modification Sites

`site.parquet` holds **protein modification sites**: one row per claim that a residue of a protein carries a modification. A site carries three separate measures of confidence:

- **`site_qvalue`**: the estimated false share among sites at least as well supported. It answers "is this residue modified at all?", estimated by the producer, for example from decoy or knockoff sites.
- **`occupancy`**: the probability that the residue carries a modification, from the spectra of the peptides that cover it.
- **`identity_probability`**: given that the residue is modified, the probability that the modification is this one. Mass spectrometry often cannot tell near-isobaric modifications apart (acetyl and trimethyl differ by 0.036 Da), so a residue can have several rows whose identity probabilities sum to 1. `alternatives` lists the competing modifications on each row.

The YAML schema definition lives at `qpx/core/data/schemas/site.yaml` in the source repository.

## Use Cases

- **Site-level PTM reporting**: report a modification at a protein position with its own error control, independent of the PSM q-value.
- **Ambiguous modifications**: keep both readings of a residue (acetyl 0.71, trimethyl 0.29) instead of a forced choice.
- **Prior annotations**: tell sites that confirm a known annotation (`prior_evidence`) from new ones.

## Schema

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `site_id` | Opaque int64 identity, derived from the identity composite (`protein_accession`, `position`, `mass_delta`, `grouped_runs`) | int64 | yes |
| `protein_accession` | Protein the position refers to | string | yes |
| `family` | Protein family (the canonical entry of an isoform, variant or chain) | string, null | no |
| `position` | 1-based residue position in the protein | int32 | yes |
| `amino_acid` | Residue at the position | string | yes |
| `modification_name` | Modification name (e.g. `Phospho`) | string | yes |
| `modification_accession` | Modification accession (e.g. `UNIMOD:21`) | string, null | no |
| `mass_delta` | Monoisotopic mass shift (Da) | float64 | yes |
| `grouped_runs` | Runs whose PSMs support the site | list\<string\> | yes |
| `pg_id` | Protein group holding the protein (`pg.pg_id`) | int64, null | no |
| `prior_evidence` | Evidence class of the claim before the search (e.g. `curated`, `catalogue`, `sample`) | string, null | no |
| `site_qvalue` | Site-level q-value | float64, null | no |
| `occupancy` | Probability that the residue is modified | float64, null | no |
| `identity_probability` | Probability of this modification given that the residue is modified | float64, null | no |
| `alternatives` | Other modifications at the residue with their identity probabilities, largest first | list\<site_alternative\>, null | no |
| `psm_count` | PSMs at the PSM cutoff that place the modification there | int32, null | no |
| `best_psm_id` | Best of those PSMs (`psm.psm_id`) | int64, null | no |
| `is_decoy` | Decoy or knockoff site | bool | yes |
| `additional_scores` | Additional scores, e.g. a localization probability | list\<score\> | no |
| `cv_params` | CV parameters | list\<cv_param\> | no |

The `site_alternative` struct:

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Modification name |
| `accession` | string, null | Modification accession |
| `mass_delta` | float64 | Monoisotopic mass shift (Da) |
| `identity_probability` | float64 | Its identity probability at the residue |

!!! note "Site, occupancy and localization"
    `site_qvalue` and `occupancy` are about the residue being modified; `identity_probability` is about which modification. Neither is a localization score: which residue of a peptide carries a modification is reported per PSM in `psm.modifications[].positions[].scores` (see [Modifications](modifications.md)).

## Python API

```python
import qpx

ds = qpx.open_dataset("PXD014414/")

# Sites on one protein
tp53 = ds.site.by_protein("P04637")

# Confident acetylation sites: site q <= 0.05, occupancy >= 0.5, identity >= 0.95
acetyl = ds.site.targets_only().by_modification("Acetyl").confident(0.05, min_occupancy=0.5, min_identity=0.95)
print(acetyl.to_df())
```

## Relationship to Other Views

- `pg_id` links a site to its protein group in [pg](pg.md).
- `best_psm_id` links to the [psm](psm.md) row with the best evidence.
- A [proteoform](proteoform.md) of the same `family` can carry the residue in a variant or chain.
