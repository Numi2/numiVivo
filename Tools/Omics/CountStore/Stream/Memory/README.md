# Count-stream temporary memory

The original reader in `41308789` retained Foundation temporaries across a long
input stream. A pool around each complete read/hash/consume iteration releases
those temporaries. The paired controls below retain all records and use the
actual product CLI, with exact integer expectations and byte-identical plans,
reports and stream hashes across implementations.

| Synthetic input | Original peak RSS | Corrected peak RSS |
| --- | ---: | ---: |
| 10 million records / 160 MB | 265.24 MB | 104.79 MB |
| 100 million records / 1.6 GB | 1,711.64 MB | 105.74 MB |

Both controls have 100,000 features and one count per feature per cell; they
use 100 and 1,000 cells respectively. The aggregate has only 100,000 entries in
both cases. The corrected peak grows by 950,272 bytes despite a tenfold larger
stream. All cell totals and aggregate counts match their analytic expectations;
all four native regression tests pass again. RSS is measured in bytes using
macOS `wait4` on each child, and the table uses decimal MB.

These are allocator/numerical controls, **not biological validation, a general
memory maximum, a throughput comparison or complete Parse qualification**.
Metadata, QC and aggregates still occupy memory. The input read is at most
1 MiB; its chunk can coexist with a parser buffer of at most 1 MiB plus 15 bytes.
The initial full Parse stream stopped on an upstream HTTP 500 after
1,166,913,217 records; the native owner rejected truncation and published no
complete bundle. Its failure is retained. The corrected owner now processes
all twelve donors as disjoint complete partitions, with transient-source retries
and independently checked native donor receipts. Full-cohort ingestion and
replay remain pending; successful donors can be resumed without reprocessing.

```sh
OPENBLAS_NUM_THREADS=1 python Tools/Omics/CountStore/Stream/Memory/run.py \
  --binary /path/to/numivivo-omics --out /path/to/control-results --label corrected
```

[Retained controls and identities](evidence/2026-09-11/manifest.json) include both
qualified executables, all original/final plans and native outputs, source hashes,
logs and the comparison. Original executable SHA-256:
`146cba39ff175a16191b0cda5a65cede65cbf88fdab18aa11f65314b388ec347`.
Corrected executable SHA-256:
`20b13e585e1526dba7be566484373a3e4775efebd66f22b3b20dff31c1ee3516`.
