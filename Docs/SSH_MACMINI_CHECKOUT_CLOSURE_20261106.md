# SSH Mac mini checkout closure — 6 November 2026

Archived the available working state for `Numi2/numiVivo`. This is a recovery checkpoint, not a qualified feature release. No simulation, build, training or experiment was resumed.

This owner has **3 registered locations**, **1 dirty checkouts** and **0 missing registrations**. The complete Mac mini audit found 123 existing Git locations, 94 dirty; there were no active project jobs or merge/rebase operations. Local source trees and staging state were preserved. Only this note and its recovery index were added to main.

No runtime features or scientific results were changed or requalified.

## Recovery

- [Release and archived payloads](https://github.com/Numi2/numiVivo/releases/tag/ssh-macmini-checkpoints-20261106).
- [Compact recovery index](SSH_MACMINI_CHECKOUT_CLOSURE_20261106.json): exact SHA-256 values, payload references, full inventory member and verification receipts.
- [Full checkpoint inventory archive](https://github.com/Numi2/numiVivo/releases/download/ssh-macmini-checkpoints-20261106/numiVivo-ssh-macmini-20261106-index.tar.zst): source-store refs/reflogs, staging state, patches and original file-to-archive mappings. Changed/untracked paths owned here: **36**.

Every new archive was independently downloaded, decompressed and restored member-by-member with size and SHA-256 checks. Git bundles were independently unpacked and all expected new objects checked in the restored packs. Two shallow third-party histories required public upstream ancestry recovery: OpenCourant for Lab and PepQuery for Vivo. Their identities and original origins are retained; no third-party repository was modified.

Fetch release tag `ssh-macmini-checkpoints-20261106` from this repository for a bundle's published-history prerequisites. Download the index archive and manifest, verify the listed SHA-256, then decompress using `zstd -d --long=29`. The inventory maps original paths to standard tar members in payload archives. Restore into a separate directory; use recorded HEADs, bundle refs, patches and index entries to reconstruct the checkpoint. Do not launch archived queues or experiments automatically.

Original heads, refs, index entries, working status, file metadata/hashes and symlink targets were rechecked before publication. Untested edits and old proposals remain unvalidated. Ignored build/cache outputs and unrelated run directories were not swept into this checkout snapshot; existing scientific evidence and the October recovery store remain untouched. There was no reset, cleanup or deletion of local source trees.
