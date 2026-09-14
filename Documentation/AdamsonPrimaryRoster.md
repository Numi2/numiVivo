# Adamson primary control roster

The primary Adamson supplement was retrieved directly on 2026-09-14 as
[Table S1](https://ars.els-cdn.com/content/image/1-s2.0-S0092867416316609-mmc1.xlsx).
The 48,003-byte workbook is `mmc1.xlsx`, SHA-256
`9b5935cb15ba2f6d60d3017832de2918e7d4f172db6f202be7999cba5feea82b`.
It contains 102 data rows: 91 rows with a `Perturb-seq_Vector_ID`, 98 rows
with a gene label (some labels cover more than one gene), and four explicit
non-targeting sequences. Its columns are `Gene`, `Protospacer`,
`Guide_ID (synonymous with sgGuide_ID)`, and `Perturb-seq_Vector_ID`.
The [retrieval manifest](../Tools/Omics/PerturbationPrediction/Adamson/evidence/2026-09-14-table-s1-retrieval/manifest.json)
and [inventory record](../Tools/Omics/PerturbationPrediction/Adamson/evidence/2026-09-14-table-s1-retrieval/table-s1.json)
retain the request log and the bounded content check; the workbook itself is
retained externally by this hash.

The paper's methods make the UPR screen's control set explicit: the large
UPR Perturb-seq experiment included **NegCtrl-2 and NegCtrl-3**. Their
protospacers are present in Table S1 and the methods identify their sequences
as follows:

| Primary label | Protospacer | Scope |
| --- | --- | --- |
| NegCtrl-2 | `GCGATGGGGGGGTGGGTAGC` | UPR Perturb-seq control |
| NegCtrl-3 | `GACGACTAGTTAGGCGTGTA` | UPR Perturb-seq control |
| NegCtrl-1 | `GGCCAAACGTGCCCTGACGG` | Other control use in the paper |
| NegCtrl-hU6 | `GCCTTGGCTAAACCGCTCCC` | Other control use in the paper |

The [paper methods](https://escholarship.org/content/qt9961h6m4/qt9961h6m4_noSplash_f242fbdc7b9f06f490340b23dde1a1ae.pdf)
also state that NegCtrl-2 was represented at eight-fold excess and that
NegCtrl-3 was used as the control in the SEC61 and HSPA5 analyses. This
confirms the experiment-level control identities and the fixed sequence
roster without using expression outcomes.

The restored GEO records still use deposited labels that are absent from
Table S1's vector column; an exhaustive value check found no `pBA580`,
`pBA582`, or `pBA581` entry:

| Restored source label | Selected cells | Current interpretation |
| --- | ---: | --- |
| `63(mod)_pBA580` | 4,595 | Candidate for one UPR control; label-to-sequence mapping not supplied |
| `Gal4-4(mod)_pBA582` | 646 | Candidate for one UPR control; label-to-sequence mapping not supplied |
| `62(mod)_pBA581` | 0 | No selected cells in the frozen UPR cohort |

These are the two selected non-gene labels in the restored cohort, but their
mapping to NegCtrl-2 versus NegCtrl-3 is not stated in the primary supplement
or the GEO guide CSV. The source also contains 94 selected guide groups while
the paper summarizes 93 guides for this experiment. The role gate therefore
remains **partially closed**: the primary control set and sequences are known,
while deposited-label mapping and the 94-versus-93 roster discrepancy still
need an author record. The separately declared sensitivity protocol below
provides conditional evidence without resolving that authoritative identity
gap.

No control assignment is inferred from abundance, the label spelling, or an
expression profile. The authoritative Adamson predictor remains unfitted and
unscored under the frozen role protocol. A separately declared
[provisional role-sensitivity experiment](AdamsonProvisionalRoleSensitivity.md)
now reports the conditional result obtained by pooling the two candidate labels;
it does not resolve their mapping or qualify the predictor.
