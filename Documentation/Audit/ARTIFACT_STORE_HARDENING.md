# Artifact store audit changes

The artifact store now uses a pinned root directory and opens every descendant directory with `openat`, `O_DIRECTORY`, `O_CLOEXEC`, and `O_NOFOLLOW`. Reads reject special files, bound allocation before reading, handle interrupted reads, and detect size changes. Immutable objects/descriptors publish with a no-clobber hard-link operation; mutable references publish with `renameat`. Temporary files use mode 0600; newly created directories use 0700. File data is fsynced before publication. Directory fsync is attempted where supported; filesystem-specific crash durability is not guaranteed by this code alone.

The root URL is a trusted caller-selected location. Its final component cannot be a symlink. The open descriptor remains authoritative if the directory is renamed. This is not a sandbox against another process with the same OS privileges, hostile mounts, or administrators. Object hashes still verify content; hashes and embedded public keys do not establish source trust without a separate trust policy.

Descriptor reads bind the requested fingerprint, canonical relative object path, length bound and metadata format. Named references must match the actual persisted descriptor, not merely supply a valid object hash alongside forged metadata. Reference reads also verify the object. Identical bytes with conflicting immutable metadata fail explicitly; encode a distinct provenance envelope when different metadata is required.

New reference filenames hash the exact UTF-8 name under `refs/v2`. This avoids case-folding and Unicode filename aliases. Valid v1 references remain readable through a checked, read-only fallback. Deletion writes a v2 tombstone so an old v1 entry cannot resurrect a deleted logical reference. Stored dates use the existing canonical-JSON second precision; put/set return the decoded persisted form rather than a subtly different subsecond in-memory value.

Default I/O bounds are 512 MiB per object, 1 MiB per descriptor/reference and 100,000 entries per listing. Callers may explicitly configure positive limits. Listing traverses only the expected SHA-256 shard hierarchy and does not silently treat malformed or inaccessible stores as complete empty listings.

## Portable check performed

The actual filesystem helper was typechecked with Swift 6 on Linux and executed against 20 targeted checks. These covered regular/empty round trips, immutable no-clobber behavior, bounded reads/listing, traversal rejection, file and directory symlinks, FIFO rejection, safe replacement of a symlink entry, and root rename/replacement while a descriptor remains open. These checks did not exercise CryptoKit, the complete artifact-store actor, APFS/iOS behavior, crash recovery or process concurrency.

Reproduce the narrow check without building NumiVivo or Metal:

```sh
scratch=$(mktemp -d)
swiftc -swift-version 6 \
  Sources/NumiVivoKit/Artifacts/VivoRootedFileStore.swift \
  Tools/Audit/RootedFileStoreChecks.swift -o "$scratch/rooted-checks"
"$scratch/rooted-checks"
rm -rf "$scratch"
```

Reference: Apple's `open(2)` documentation describes O_NOFOLLOW, O_EXCL and nonblocking opens: https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man2/open.2.html
