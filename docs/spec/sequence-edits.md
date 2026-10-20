# Sequence Edits

A peptide from a sequence variant differs from the reference protein: a single substitution, several adjacent residues replaced, residues deleted or inserted. The [peptidoform](peptidoform.md) records only the observed sequence, and ProForma has no notation for a change relative to a reference, so the change itself would otherwise not be recorded anywhere. QPX records it as a **sequence edit**: a span of the reference protein and what replaces it.

## The `sequence_edit` struct

```text
sequence_edit: struct{
    start:       int32,          -- first reference residue of the span (1-based)
    end:         int32,          -- last reference residue of the span (inclusive)
    reference:   string,         -- the reference residues of start..end ("" for an insertion)
    alternative: string,         -- what replaces them ("" for a deletion)
    name:        string, null,   -- e.g. HGVS protein notation "p.Arg175His"
    accession:   string, null    -- e.g. a dbSNP, COSMIC or UniProt variant identifier
}
```

Coordinates are the reference protein's, 1-based and inclusive, so `reference` holds exactly `end - start + 1` residues. An insertion has an empty `reference` and `end = start - 1`: it goes before residue `start`.

| Change | start | end | reference | alternative | name |
| ------ | ----- | --- | --------- | ----------- | ---- |
| Substitution | 175 | 175 | `R` | `H` | p.Arg175His |
| Two adjacent residues replaced | 10 | 11 | `KE` | `QK` | p.Lys10_Glu11delinsGlnLys |
| Deletion | 30 | 30 | `G` | `` | p.Gly30del |
| Insertion before residue 50 | 50 | 49 | `` | `AA` | p.Pro49_Leu50insAlaAla |

Several edits on one molecule are several entries, sorted by `start` and not overlapping.

## Where edits are stored

| View | Field | Struct | Meaning |
| ---- | ----- | ------ | ------- |
| [pepmap](pepmap.md) | `sequence_edits` | `peptide_sequence_edit` (a `sequence_edit` with `protein_accession`) | The changes a peptide carries relative to a reference protein it maps to; null when it matches its proteins as they are |
| [proteoform](proteoform.md) | `edits` | `sequence_edit` | For a variant member, the changes that make it from the canonical entry, in canonical coordinates |

Both columns are optional.

## Substitutions and modifications

Some substitutions have the same mass as a modification: D→E and G→A weigh the same as methylation, and N→D the same as deamidation. A spectrum that fits both readings is a case of [competing interpretations](modifications.md#competing-interpretations): a single-residue substitution can be listed as a candidate at its position, using Unimod's amino-acid substitution entries (e.g. `Asp->Glu`), next to the modification it competes with. Edits that span several residues cannot sit at one position and are recorded only as sequence edits.

## Validation

Under `qpxc validate` (strict), each edit has `start >= 1`, a `reference` of `end - start + 1` one-letter residues (empty with `end = start - 1` for an insertion), an `alternative` of one-letter residues, and changes something (not both empty).
