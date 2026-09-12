# HIRISA donor-level exposure join

All five frozen HIRISA training donors join exactly to their B-cell IFNb and
unstimulated libraries in the retained GSE306664 design: ten libraries total.
The executable check verifies donor, cell population, treatment, duration, batch
and pool for both members of each pair. Each retained sample protocol specifies
IFN-beta at 100 units/mL, R&D catalogue 8499-IF-010/CF, and a 21-hour incubation.
The control condition is explicitly `none`; it is not a fresh-cell control.

`mapping.json` retains the five model IDs, ten accessions, selected source fields
and original design/model metadata hashes, plus the upstream SOFT hash recorded
by the design. `map.py` performs the join against those retained files. This turn
did not download the original SOFT again or reprocess expression values. A live
[primary sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9205455)
also confirms the protocol; four other live sample fetches returned cache misses,
so their verification uses the existing deposited-metadata artifact.

This resolves HIRISA donor/library metadata linkage and the vendor/catalogue gap
left by the shorter institute methods page. It does not convert activity units to
mass concentration, establish reagent-lot equivalence, participant independence
across studies, or identify a causal duration effect. The matching catalogue string
in Parse is not proof of matched exposure. Model fitting and outcomes are unchanged.
