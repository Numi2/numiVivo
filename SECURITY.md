# Security policy

NumiVivo processes untrusted program documents, biological model files, mechanism libraries, and experiment artifacts. Treat all imported content as hostile until it passes structural, semantic, and resource-bound validation.

## Supported versions

Security fixes are applied to the current `main` branch until tagged releases are introduced. Artifact readers must reject unsupported major versions rather than attempting permissive recovery.

## Report privately

Do not publish an exploit, unsafe generated construct, or patient-associated data in a public issue. Send a private GitHub security advisory to the repository maintainers with:

- affected revision;
- input artifact or minimal reproducer;
- observed impact;
- whether the issue crosses the computational/physical boundary;
- suggested containment, when known.

## Security invariants

The following are release-blocking defects:

- out-of-bounds section, table, string, or GPU-buffer access;
- integer overflow in counts, offsets, dispatch dimensions, or allocation sizes;
- unbounded parser recursion or attacker-controlled allocation;
- acceptance of overlapping or unauthenticated pack sections;
- state publication after a failed numerical or safety check;
- nondeterministic reuse of random streams across cells or experiments;
- silent unit conversion, identifier substitution, or mechanism fallback;
- loss of provenance linking output to exact input fingerprints;
- execution of embedded scripts, shell commands, network requests, or dynamic libraries from an artifact;
- exposure of personally identifying, genomic, clinical, or regulated data through diagnostics or telemetry.

## Computational safety boundary

NumiVivo does not authorize clinical, environmental, or in-vivo use. A generated molecular design or favourable simulation result is not a safety determination. Physical implementation requires independent review, containment, ethics approval, and applicable regulatory authorization.

The repository does not accept secrets, patient records, identifiable genomic data, access tokens, or proprietary sequence libraries in examples, issues, logs, or fixtures.
