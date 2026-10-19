# Proteoforms

`proteoform.parquet` reports, for each protein **family**, the canonical entry and each **isoform**, **sequence variant** and **processed chain** that the data support with their own evidence. A family groups a canonical protein (e.g. a UniProt entry) with its alternative forms, so alternative entries do not split the protein's evidence; each member is a separate claim with its own error control.

The YAML schema definition lives at `qpx/core/data/schemas/proteoform.yaml` in the source repository.

## Use Cases

- **Proteogenomics**: report which sequence variants of a protein are observed, with a member-level q-value.
- **Isoform and chain evidence**: separate the evidence specific to an isoform or a processed chain from the evidence shared with the canonical form.
- **Proteoform quantification**: quantify a member from its specific peptides, per run or channel.

## Schema

| Field | Description | Type | Required |
|-------|-------------|------|----------|
| `proteoform_id` | Opaque int64 identity, derived from the identity composite (`accession`, `grouped_runs`, `label`) | int64 | yes |
| `family` | Canonical entry shared by all members of the family | string | yes |
| `accession` | This member's accession | string | yes |
| `kind` | `canonical`, `isoform`, `variant` or `chain` | string | yes |
| `description` | What distinguishes the member (e.g. `p.Arg175His`, an isoform identifier, a chain range) | string, null | no |
| `pg_id` | Protein group holding the family (`pg.pg_id`) | int64, null | no |
| `grouped_runs` | Runs grouped into this quantification unit | list\<string\> | yes |
| `source` | Database or annotation the member comes from | string, null | no |
| `prior_evidence` | Evidence class of the member before the search (e.g. `curated`, `catalogue`, `sample`) | string, null | no |
| `specific_peptides` | Peptides no other member of the family explains | int32, null | no |
| `shared_peptides` | Peptides shared with other members of the family | int32, null | no |
| `psm_count` | PSMs at their cutoff that support the member | int32, null | no |
| `evidence` | `specific` or `shared`: where `qvalue` comes from | string, null | no |
| `qvalue` | Member-level q-value | float64, null | no |
| `best_psm_id` | PSM the q-value comes from (`psm.psm_id`) | int64, null | no |
| `is_decoy` | Decoy or knockoff member | bool | yes |
| `label` | Channel/label of the quantity; null for identification only | string, null | no |
| `intensity` | Quantity from the member's specific peptides for this label | float32, null | no |
| `additional_intensities` | Other quantities for this label | list\<additional_intensity\> | no |
| `additional_scores` | Additional scores | list\<score\> | no |
| `cv_params` | CV parameters | list\<cv_param\> | no |

!!! note "One row per label"
    Like [pg](pg.md), the view is flattened: a quantified member has one row per label, and an identification-only member has one row with `label` and `intensity` null. Writers accept the natural `intensities: list<{label, intensity}>` per member and unroll it.

## Python API

```python
import qpx

ds = qpx.open_dataset("PXD014414/")

# Variants observed at a member q-value of 0.05
variants = ds.proteoform.targets_only().by_kind("variant").confident(0.05)
print(variants.to_df()[["family", "accession", "description", "qvalue"]])

# All members of one family
print(ds.proteoform.by_family("P04637").to_df())
```

## Relationship to Other Views

- `pg_id` links the family to its protein group in [pg](pg.md); the pg q-value controls the protein, the proteoform q-value the member's own claim.
- `best_psm_id` links to the [psm](psm.md) row the q-value comes from.
- [site](site.md) rows of the same `family` report modification sites.
