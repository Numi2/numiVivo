# Eligible-level index follow-up

Declared 2026-09-11 after the first candidate completed Hagai/Kang fitting and
replay but failed Ding's fixed 500-million-evaluation index work limit. No
biological metrics have been read. Preserve the first candidate, its source,
protocol, completed outputs, incomplete Ding output and failed status.

Keep every numerical/biological setting and evaluation rule in PROTOCOL.md,
except replace the single globally filtered matching index as follows. Build
one incremental lower-level index in ascending level order, and one incremental
upper-level index in descending level order. Within each level, insert original
source rows in ascending order. Query a level before inserting its cells; thus
only eligible other levels are in its index. Skip insertion of the final level
in each direction, whose cells cannot be references for any remaining query.
Each index uses M=16, construction ef=200, search ef=128, seed=7 and the unchanged
500-million-metric-evaluation lifetime budget. Process the two indexes serially;
report their combined construction/query work and maximum serialized index size.
This changes approximate search topology and is a new candidate, not a resume
or a relabeling of first-attempt outputs.

The 64-anchor local kernel, all source cells and anchor records, sigma, scaling,
panorama rules, seed, fixed biological margins, complete replay and independent
checks remain unchanged. Refit all three original cohorts; do not retain only
the previously successful cohorts. This is one design repair to filtered-search
work admission, not a hyperparameter or outcome-based search.
