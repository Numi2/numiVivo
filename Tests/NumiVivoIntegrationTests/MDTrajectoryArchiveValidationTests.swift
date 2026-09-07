import Foundation
import Darwin
import Testing
@testable import NumiVivoKit

/// Storage-only checks: none of these tests constructs a native MD runtime.
@Suite(.serialized) struct MDTrajectoryArchiveValidationTests: Sendable {
    private struct Fixture: Sendable {
        let root: URL
        let store: VivoArtifactStore
        let manifest: VivoStoredArtifact
        let reader: VivoMDTrajectoryArchiveReader
        let links: [VivoStoredArtifact]
        let payloads: [VivoStoredArtifact]
        let frames: [[VivoMDArchiveFrame]]
    }

    private func frame(_ ordinal: UInt64, velocities: Bool) -> VivoMDArchiveFrame {
        .init(stepIndex: 10 + ordinal * 2, timePS: 2 + Double(ordinal) * 0.25,
              positionsNM: [.init(Double(ordinal % 16) * 0.125, 0, -0.5)],
              velocitiesNMPerPS: velocities ? [.init(0.25, -0.5, 1)] : nil, periodicCell: nil)
    }

    private func storeJSON<T: Encodable & Sendable>(_ value: T, kind: String,
                                                   store: VivoArtifactStore) async throws -> VivoStoredArtifact {
        try await store.put(data: VivoCanonicalJSON.encode(value), kind: kind,
                            mediaType: "application/vnd.numivivo.\(kind)+json")
    }

    private func fixture(_ chunkFrames: [Int] = [2, 1, 2], velocities: Bool = true) async throws -> Fixture {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-archive-validation-\(UUID().uuidString)")
        do {
            let store = try VivoArtifactStore(rootURL: root)
            let system = try VivoCanonicalJSON.fingerprint(Data("trajectory validation system".utf8))
            let configuration = try VivoCanonicalJSON.fingerprint(Data("trajectory validation configuration".utf8))
            var ordinal: UInt64 = 0
            var links: [VivoStoredArtifact] = [], payloads: [VivoStoredArtifact] = []
            var chunks: [[VivoMDArchiveFrame]] = []
            for count in chunkFrames {
                let frames = (0..<count).map { frame(ordinal + UInt64($0), velocities: velocities) }
                let first = try #require(frames.first), last = try #require(frames.last)
                let bytes = try VivoMDTrajectoryChunkCodec.encode(frames, systemFingerprint: system,
                    configurationFingerprint: configuration, particleCount: 1,
                    includeVelocities: velocities, firstFrameOrdinal: ordinal)
                let payload = try await store.put(data: bytes, kind: "md-trajectory-chunk",
                                                  mediaType: VivoMDTrajectoryChunkCodec.mediaType)
                let link = VivoMDTrajectoryChunkLink(schema: VivoMDTrajectoryChunkLink.schemaID,
                    systemFingerprint: system, configurationFingerprint: configuration,
                    particleCount: 1, includeVelocities: velocities,
                    chunkOrdinal: UInt64(links.count), firstFrameOrdinal: ordinal,
                    frameCount: UInt32(count), previous: links.last?.fingerprint,
                    payload: payload.fingerprint, payloadBytes: payload.byteCount,
                    firstStep: first.stepIndex, lastStep: last.stepIndex,
                    firstTimePS: first.timePS, lastTimePS: last.timePS)
                links.append(try await storeJSON(link, kind: "md-trajectory-link", store: store))
                payloads.append(payload); chunks.append(frames); ordinal += UInt64(count)
            }
            let first = chunks.first?.first, last = chunks.last?.last
            let value = VivoMDTrajectoryManifest(schema: VivoMDTrajectoryManifest.schemaID,
                systemFingerprint: system, configurationFingerprint: configuration,
                particleCount: 1, includeVelocities: velocities, frameCount: ordinal,
                chunkCount: UInt64(links.count), tail: links.last?.fingerprint,
                firstStep: first?.stepIndex, lastStep: last?.stepIndex,
                firstTimePS: first?.timePS, lastTimePS: last?.timePS, sealed: false)
            let manifest = try await storeJSON(value, kind: "md-trajectory-manifest", store: store)
            let reader = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: manifest.fingerprint)
            return .init(root: root, store: store, manifest: manifest, reader: reader,
                         links: links, payloads: payloads, frames: chunks)
        } catch {
            try? FileManager.default.removeItem(at: root)
            throw error
        }
    }

    private func altered<T: Codable>(_ value: T, _ edit: (inout [String: Any]) -> Void) throws -> T {
        var object = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(value)) as? [String: Any])
        edit(&object)
        return try VivoCanonicalJSON.decode(T.self, from: JSONSerialization.data(withJSONObject: object, options: [.sortedKeys]))
    }

    /// Rehash every successor so structural tests reach the changed contract,
    /// rather than failing merely because its content-addressed bytes changed.
    private func rewritten(_ f: Fixture, edit: (Int, inout [String: Any]) -> Void) async throws -> VivoMDTrajectoryArchiveReader {
        var previous: VivoFingerprint?
        for (ordinal, record) in f.links.enumerated() {
            let original = try await f.reader.readLink(record.fingerprint)
            let value = try altered(original) { object in
                if let previous { object["previous"] = ["bytes": previous.bytes] }
                else { object.removeValue(forKey: "previous") }
                edit(ordinal, &object)
            }
            previous = try await storeJSON(value, kind: "md-trajectory-link", store: f.store).fingerprint
        }
        let manifest = try altered(f.reader.manifest) { object in
            if let previous { object["tail"] = ["bytes": previous.bytes] }
        }
        let stored = try await storeJSON(manifest, kind: "md-trajectory-manifest", store: f.store)
        return try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: stored.fingerprint)
    }

    private func expectInvalid(_ operation: () async throws -> Void) async {
        do { try await operation(); Issue.record("Expected an invalid trajectory contract") }
        catch is VivoArtifactValidationError { }
        catch { Issue.record("Expected trajectory validation failure, received \(error)") }
    }

    private func expectReadBound(_ operation: () async throws -> Void) async {
        do { try await operation(); Issue.record("Expected the rooted file read to reject its preallocation bound") }
        catch VivoRootedFileStore.Failure.exceededLimit(_) { }
        catch { Issue.record("Expected rooted read bound failure, received \(error)") }
    }

    private func expectCancellation(_ operation: () async throws -> Void) async {
        do { try await operation(); Issue.record("Cancellation returned a successful partial result") }
        catch is CancellationError { }
        catch { Issue.record("Expected CancellationError, received \(error)") }
    }

    private func expectIntegrityFailure(_ fingerprint: VivoFingerprint,
                                        _ operation: () async throws -> Void) async {
        do { try await operation(); Issue.record("Expected corrupt content to fail integrity verification") }
        catch VivoArtifactStoreError.integrityFailure(let actual) { #expect(actual == fingerprint) }
        catch { Issue.record("Expected payload integrity failure, received \(error)") }
    }

    @Test(arguments: [[], [1], [1, 1, 1], [2, 1, 3]])
    func streamingScopesMatchChronologicalFrames(chunkFrames: [Int]) async throws {
        for velocities in [false, true] {
            let f = try await fixture(chunkFrames, velocities: velocities)
            defer { try? FileManager.default.removeItem(at: f.root) }
            #expect(try await f.reader.index() == f.links.map(\.fingerprint))
            var decoded: [VivoMDArchiveFrame] = []
            for link in try await f.reader.index() { decoded += try await f.reader.readChunk(link) }
            #expect(decoded == f.frames.flatMap { $0 })
            for scope in [VivoMDTrajectoryValidationScope.index, .restart, .allPayloads] {
                let result = try await f.reader.validate(scope: scope)
                let expectedPayloads = scope == .index ? [] : scope == .restart ? Array(f.payloads.suffix(1)) : f.payloads
                #expect(result.manifestFingerprint == f.manifest.fingerprint)
                #expect(result.scope == scope && result.indexedChunks == UInt64(chunkFrames.count))
                #expect(result.verifiedPayloads == UInt64(expectedPayloads.count))
                #expect(result.verifiedPayloadBytes == expectedPayloads.reduce(UInt64(0)) { $0 + $1.byteCount })
            }
        }
    }

    @Test func explicitCallerLimitsRemainDistinctFromUncappedValidation() async throws {
        let f = try await fixture([1, 1, 1, 1, 1])
        defer { try? FileManager.default.removeItem(at: f.root) }
        #expect(try await f.reader.validate(scope: .allPayloads).indexedChunks == 5)
        #expect(try await f.reader.index(maximumChunks: 5) == f.links.map(\.fingerprint))
        for limit in [UInt64(0), 4] {
            await expectInvalid { _ = try await f.reader.index(maximumChunks: limit) }
            for scope in [VivoMDTrajectoryValidationScope.index, .restart, .allPayloads] {
                await expectInvalid { _ = try await f.reader.validate(scope: scope, maximumChunks: limit) }
            }
        }
        #expect(try await f.reader.validate(scope: .restart, maximumChunks: 5).verifiedPayloads == 1)
        // A caller work limit is not evidence that a manifest's declared count
        // is real. The first link must reject this forged ordinal immediately,
        // before index() tries to allocate UInt64(Int.max) array elements.
        let forged = try altered(f.reader.manifest) { object in
            object["chunkCount"] = UInt64(Int.max)
            object["frameCount"] = UInt64(Int.max)
        }
        let stored = try await storeJSON(forged, kind: "md-trajectory-manifest", store: f.store)
        let reader = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: stored.fingerprint)
        do {
            _ = try await reader.index(maximumChunks: UInt64(Int.max))
            Issue.record("A forged enormous manifest count was accepted")
        } catch VivoArtifactValidationError.invalid(let message) {
            #expect(message == "trajectory link ordinal gap")
        }
    }

    @Test func olderPayloadCorruptionIsOutsideRestartScopeButInsideExhaustiveVerification() async throws {
        let f = try await fixture()
        defer { try? FileManager.default.removeItem(at: f.root) }
        let reference = try await f.store.setReference("archive-prefix", to: f.manifest)
        let old = try #require(f.payloads.first)
        var bytes = try await f.store.data(for: old.fingerprint)
        bytes[bytes.startIndex] ^= 0xff
        try bytes.write(to: f.root.appendingPathComponent(old.objectPath))
        #expect(try await f.reader.validate(scope: .index).verifiedPayloads == 0)
        let restart = try await f.reader.validate(scope: .restart)
        #expect(restart.indexedChunks == 3 && restart.verifiedPayloads == 1)
        #expect(restart.verifiedPayloadBytes == f.payloads.last?.byteCount)
        _ = try await VivoMDTrajectoryArchiveWriter.resume(store: f.store, manifest: f.manifest.fingerprint)
        await expectIntegrityFailure(old.fingerprint) { _ = try await f.reader.validate(scope: .allPayloads) }
        await expectIntegrityFailure(old.fingerprint) { try await f.reader.verify() }
        #expect(try await f.store.reference("archive-prefix") == reference)
    }

    @Test func corruptFinalPayloadPreventsRestartAndDoesNotChangeTheReference() async throws {
        let f = try await fixture()
        defer { try? FileManager.default.removeItem(at: f.root) }
        let reference = try await f.store.setReference("archive-prefix", to: f.manifest)
        let tail = try #require(f.payloads.last)
        var bytes = try await f.store.data(for: tail.fingerprint)
        bytes[bytes.startIndex] ^= 0xff
        try bytes.write(to: f.root.appendingPathComponent(tail.objectPath))
        #expect(try await f.reader.validate(scope: .index).indexedChunks == 3)
        for scope in [VivoMDTrajectoryValidationScope.restart, .allPayloads] {
            await expectIntegrityFailure(tail.fingerprint) { _ = try await f.reader.validate(scope: scope) }
        }
        await expectIntegrityFailure(tail.fingerprint) {
            _ = try await VivoMDTrajectoryArchiveWriter.resume(store: f.store, manifest: f.manifest.fingerprint)
        }
        #expect(try await f.store.reference("archive-prefix") == reference)
    }

    @Test(arguments: ["truncated", "magic", "ordinal", "range"])
    func rehashedMalformedPayloadStillFailsTheRestartContract(fault: String) async throws {
        let f = try await fixture()
        defer { try? FileManager.default.removeItem(at: f.root) }
        let tail = try #require(f.payloads.last)
        var bytes = try await f.store.data(for: tail.fingerprint)
        if fault == "truncated" { bytes.removeLast() }
        else if fault == "magic" { bytes[bytes.startIndex] = 0 }
        else {
            var frames = try #require(f.frames.last)
            if fault == "range" {
                let last = try #require(frames.last)
                frames[frames.count - 1] = .init(stepIndex: last.stepIndex, timePS: last.timePS + 0.125,
                    positionsNM: last.positionsNM, velocitiesNMPerPS: last.velocitiesNMPerPS,
                    periodicCell: last.periodicCell)
            }
            bytes = try VivoMDTrajectoryChunkCodec.encode(frames,
                systemFingerprint: f.reader.manifest.systemFingerprint,
                configurationFingerprint: f.reader.manifest.configurationFingerprint,
                particleCount: 1, includeVelocities: true, firstFrameOrdinal: fault == "ordinal" ? 999 : 3)
        }
        let replacement = try await f.store.put(data: bytes, kind: "md-trajectory-chunk",
                                                mediaType: VivoMDTrajectoryChunkCodec.mediaType)
        let reader = try await rewritten(f) { ordinal, object in
            if ordinal == 2 { object["payload"] = ["bytes": replacement.fingerprint.bytes] }
        }
        #expect(try await reader.validate(scope: .index).indexedChunks == 3)
        for scope in [VivoMDTrajectoryValidationScope.restart, .allPayloads] {
            await expectInvalid { _ = try await reader.validate(scope: scope) }
        }
    }

    private enum LinkFault: String, CaseIterable, Sendable {
        case schema, systemIdentity, configurationIdentity, particleCount, velocities
        case chunkOrdinalGap, frameOrdinalGap, frameOrdinalOverflow, zeroFrames, oversizedFrameCount, wrongPayloadBytes
        case stepOverlap, timeOverlap, wrongFirstBoundary, wrongLastBoundary, truncatedChain
    }

    @Test(arguments: LinkFault.allCases)
    private func rehashedStructuralCorruptionFailsEveryScope(fault: LinkFault) async throws {
        let f = try await fixture()
        defer { try? FileManager.default.removeItem(at: f.root) }
        let reader = try await rewritten(f) { ordinal, object in
            switch fault {
            case .wrongFirstBoundary where ordinal == 0: object["firstStep"] = 11
            case .wrongLastBoundary where ordinal == 2: object["lastTimePS"] = 2.875
            case .truncatedChain where ordinal == 2: object.removeValue(forKey: "previous")
            case .schema where ordinal == 1: object["schema"] = "unknown"
            case .systemIdentity where ordinal == 1: object["systemFingerprint"] = ["bytes": Array(repeating: UInt8(0x5a), count: 32)]
            case .configurationIdentity where ordinal == 1: object["configurationFingerprint"] = ["bytes": Array(repeating: UInt8(0x5a), count: 32)]
            case .particleCount where ordinal == 1: object["particleCount"] = 2
            case .velocities where ordinal == 1: object["includeVelocities"] = false
            case .chunkOrdinalGap where ordinal == 1: object["chunkOrdinal"] = 2
            case .frameOrdinalGap where ordinal == 1: object["firstFrameOrdinal"] = 1
            case .frameOrdinalOverflow where ordinal == 1: object["firstFrameOrdinal"] = UInt64.max
            case .zeroFrames where ordinal == 1: object["frameCount"] = 0
            case .oversizedFrameCount where ordinal == 1: object["frameCount"] = 4097
            case .wrongPayloadBytes where ordinal == 1: object["payloadBytes"] = 1
            case .stepOverlap where ordinal == 1: object["lastStep"] = 16
            case .timeOverlap where ordinal == 1: object["lastTimePS"] = 2.75
            default: break
            }
        }
        for scope in [VivoMDTrajectoryValidationScope.index, .restart, .allPayloads] {
            await expectInvalid { _ = try await reader.validate(scope: scope) }
        }
        await expectInvalid { _ = try await reader.index() }
    }

    @Test func missingLinkFailsInsteadOfReturningAPartialIndex() async throws {
        let f = try await fixture()
        defer { try? FileManager.default.removeItem(at: f.root) }
        let missing = f.links[1]
        try FileManager.default.removeItem(at: f.root.appendingPathComponent(missing.objectPath))
        for scope in [VivoMDTrajectoryValidationScope.index, .restart, .allPayloads] {
            do {
                _ = try await f.reader.validate(scope: scope)
                Issue.record("A missing predecessor was accepted")
            } catch VivoArtifactStoreError.objectMissing(let fingerprint) { #expect(fingerprint == missing.fingerprint) }
            catch { Issue.record("Expected missing-link failure, received \(error)") }
        }
    }

    @Test func rehashedInvalidManifestIsRejectedBeforeTraversal() async throws {
        let f = try await fixture()
        defer { try? FileManager.default.removeItem(at: f.root) }
        for fault in ["schema", "counts", "tail", "range", "emptyMetadata"] {
            let value = try altered(f.reader.manifest) { object in
                switch fault {
                case "schema": object["schema"] = "unknown"
                case "counts": object["chunkCount"] = 6
                case "tail": object.removeValue(forKey: "tail")
                case "range": object["firstStep"] = 100
                default: object["frameCount"] = 0; object["chunkCount"] = 0
                }
            }
            let record = try await storeJSON(value, kind: "md-trajectory-manifest", store: f.store)
            await expectInvalid { _ = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: record.fingerprint) }
        }
    }

    @Test func metadataAndPayloadBoundsAreEnforcedByTheRootedRead() async throws {
        let f = try await fixture()
        defer { try? FileManager.default.removeItem(at: f.root) }
        var manifestBytes = try await f.store.data(for: f.manifest.fingerprint)
        manifestBytes.append(Data(repeating: 0x20, count: 64 * 1024 + 1 - manifestBytes.count))
        let oversizedManifest = try await f.store.put(data: manifestBytes, kind: "md-trajectory-manifest",
                                                      mediaType: "application/vnd.numivivo.md-trajectory-manifest+json")
        await expectReadBound {
            _ = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: oversizedManifest.fingerprint)
        }

        let originalTail = try #require(f.links.last)
        var linkBytes = try await f.store.data(for: originalTail.fingerprint)
        linkBytes.append(Data(repeating: 0x20, count: 64 * 1024 + 1 - linkBytes.count))
        let oversizedLink = try await f.store.put(data: linkBytes, kind: "md-trajectory-link",
                                                  mediaType: "application/vnd.numivivo.md-trajectory-link+json")
        await expectReadBound { _ = try await f.reader.readLink(oversizedLink.fingerprint) }

        let payload = try #require(f.payloads.last)
        var payloadBytes = try await f.store.data(for: payload.fingerprint)
        payloadBytes.append(0)
        let oversizedPayload = try await f.store.put(data: payloadBytes, kind: "md-trajectory-chunk",
                                                     mediaType: VivoMDTrajectoryChunkCodec.mediaType)
        let reader = try await rewritten(f) { ordinal, object in
            if ordinal == 2 { object["payload"] = ["bytes": oversizedPayload.fingerprint.bytes] }
        }
        #expect(try await reader.validate(scope: .index).indexedChunks == 3)
        await expectReadBound { _ = try await reader.validate(scope: .restart) }
        await expectReadBound { _ = try await reader.validate(scope: .allPayloads) }
    }

    @Test func boundedArtifactReadsApplyWithVerificationDisabledAndAllowEmptyObjects() async throws {
        let f = try await fixture([])
        defer { try? FileManager.default.removeItem(at: f.root) }
        let bytes = Data(repeating: 7, count: 4096)
        let record = try await f.store.put(data: bytes, kind: "test", mediaType: "application/octet-stream")
        for verify in [false, true] {
            #expect(try await f.store.data(for: record.fingerprint, maximumBytes: bytes.count, verify: verify) == bytes)
            await expectReadBound {
                _ = try await f.store.data(for: record.fingerprint, maximumBytes: bytes.count - 1, verify: verify)
            }
        }
        let tighterStore = try VivoArtifactStore(rootURL: f.root, createIfNeeded: false,
                                                limits: .init(maximumObjectBytes: 128))
        await expectReadBound {
            _ = try await tighterStore.data(for: record.fingerprint, maximumBytes: bytes.count, verify: false)
        }
        let empty = try await f.store.put(data: Data(), kind: "empty", mediaType: "application/octet-stream")
        #expect(try await f.store.data(for: empty.fingerprint, maximumBytes: 0).isEmpty)
        do {
            _ = try await f.store.data(for: empty.fingerprint, maximumBytes: -1)
            Issue.record("A negative artifact read limit was accepted")
        } catch VivoArtifactStoreError.invalidDescriptor(_) { }
    }

    @Test func resumeTokensRequirePayloadEvidenceAndRetainTheirExactStoreAuthority() async throws {
        let first = try await fixture(), other = try await fixture()
        defer {
            try? FileManager.default.removeItem(at: first.root)
            try? FileManager.default.removeItem(at: other.root)
        }
        // Both roots contain the identical content-addressed prefix. A token
        // must still carry the authority of the actual reader that produced it.
        #expect(first.manifest.fingerprint == other.manifest.fingerprint)
        let indexOnly = try await first.reader.validate(scope: .index)
        await expectInvalid { _ = try VivoMDTrajectoryArchiveWriter.resume(validated: indexOnly) }
        let sealed = try altered(first.reader.manifest) { $0["sealed"] = true }
        let sealedArtifact = try await storeJSON(sealed, kind: "md-trajectory-manifest", store: first.store)
        let sealedReader = try await VivoMDTrajectoryArchiveReader.open(store: first.store, manifest: sealedArtifact.fingerprint)
        let sealedToken = try await sealedReader.validate(scope: .restart)
        await expectInvalid { _ = try VivoMDTrajectoryArchiveWriter.resume(validated: sealedToken) }
        for scope in [VivoMDTrajectoryValidationScope.restart, .allPayloads] {
            let token = try await first.reader.validate(scope: scope)
            let writer = try VivoMDTrajectoryArchiveWriter.resume(validated: token)
            let next = frame(first.reader.manifest.frameCount, velocities: true)
            try await writer.append(.init(systemFingerprint: first.reader.manifest.systemFingerprint,
                configurationFingerprint: first.reader.manifest.configurationFingerprint,
                stepIndex: next.stepIndex, timePS: next.timePS, positionsNM: next.positionsNM,
                velocitiesNMPerPS: try #require(next.velocitiesNMPerPS), periodicCell: nil))
            let extended = try await writer.snapshot()
            #expect(await first.store.contains(extended.fingerprint))
            #expect(await other.store.contains(extended.fingerprint) == false)
            let resumedReader = try await VivoMDTrajectoryArchiveReader.open(store: first.store, manifest: extended.fingerprint)
            let tail = try #require(resumedReader.manifest.tail)
            #expect(try await resumedReader.readLink(tail).previous == first.reader.manifest.tail)
            #expect(try await resumedReader.readChunk(tail) == [next])
        }
        let token = try await first.reader.validate(scope: .restart)
        let cancelled = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try VivoMDTrajectoryArchiveWriter.resume(validated: token)
        }
        await expectCancellation { _ = try await cancelled.value }
    }

    @Test func resumeTokenKeepsThePinnedRootWhenItsOriginalPathIsReplaced() async throws {
        let f = try await fixture()
        let movedRoot = f.root.deletingLastPathComponent().appendingPathComponent(
            "\(f.root.lastPathComponent)-moved-\(UUID().uuidString)")
        defer {
            try? FileManager.default.removeItem(at: f.root)
            try? FileManager.default.removeItem(at: movedRoot)
        }
        let token = try await f.reader.validate(scope: .restart)
        try FileManager.default.moveItem(at: f.root, to: movedRoot)
        let replacement = try VivoArtifactStore(rootURL: f.root)
        let originalRootURL = await f.store.rootURL
        let replacementRootURL = await replacement.rootURL
        #expect(originalRootURL == replacementRootURL)
        let sentinel = try await replacement.put(data: Data("replacement root".utf8), kind: "test",
                                                 mediaType: "application/octet-stream")
        #expect(await f.store.contains(sentinel.fingerprint) == false)

        // The path now names an unrelated store. Reopening token.store.rootURL
        // would write there, even though the token validated the moved inode.
        let writer = try VivoMDTrajectoryArchiveWriter.resume(validated: token)
        let next = frame(f.reader.manifest.frameCount, velocities: true)
        try await writer.append(.init(systemFingerprint: f.reader.manifest.systemFingerprint,
            configurationFingerprint: f.reader.manifest.configurationFingerprint,
            stepIndex: next.stepIndex, timePS: next.timePS, positionsNM: next.positionsNM,
            velocitiesNMPerPS: try #require(next.velocitiesNMPerPS), periodicCell: nil))
        let extended = try await writer.snapshot()
        #expect(await f.store.contains(extended.fingerprint))
        #expect(await replacement.contains(extended.fingerprint) == false)
        let movedBytes = try Data(contentsOf: movedRoot.appendingPathComponent(extended.objectPath))
        #expect(try VivoCanonicalJSON.fingerprint(movedBytes) == extended.fingerprint)
        #expect(FileManager.default.fileExists(atPath: f.root.appendingPathComponent(extended.objectPath).path) == false)
        let reader = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: extended.fingerprint)
        let tail = try #require(reader.manifest.tail)
        #expect(try await reader.readLink(tail).previous == f.reader.manifest.tail)
        #expect(try await reader.readChunk(tail) == [next])
        #expect(try await reader.validate(scope: .allPayloads).indexedChunks == 4)
        #expect(try await replacement.data(for: sentinel.fingerprint) == Data("replacement root".utf8))
    }

    @Test(arguments: [[], [1, 1, 1]])
    func alreadyCancelledValidationCannotReturnSuccess(chunkFrames: [Int]) async throws {
        let f = try await fixture(chunkFrames)
        defer { try? FileManager.default.removeItem(at: f.root) }
        let reference = try await f.store.setReference("archive-prefix", to: f.manifest)
        for scope in [VivoMDTrajectoryValidationScope.index, .restart, .allPayloads] {
            let task = Task {
                withUnsafeCurrentTask { $0?.cancel() }
                return try await f.reader.validate(scope: scope)
            }
            await expectCancellation { _ = try await task.value }
        }
        let index = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try await f.reader.index()
        }
        await expectCancellation { _ = try await index.value }
        let exhaustive = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            try await f.reader.verify()
        }
        await expectCancellation { try await exhaustive.value }
        let resume = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try await VivoMDTrajectoryArchiveWriter.resume(store: f.store, manifest: f.manifest.fingerprint)
        }
        await expectCancellation { _ = try await resume.value }
        #expect(try await f.store.reference("archive-prefix") == reference)
    }

    @Test(arguments: [false, true])
    func cancellationAfterARealTraversalReadRejectsTheNextPull(readPayload: Bool) async throws {
        let f = try await fixture()
        defer { try? FileManager.default.removeItem(at: f.root) }
        let reference = try await f.store.setReference("archive-prefix", to: f.manifest)
        let task = Task {
            var traversal = try VivoMDTrajectoryIndexTraversal(reader: f.reader)
            let first = try #require(try await traversal.next())
            #expect(first.fingerprint == f.links.last?.fingerprint)
            if readPayload { #expect(try await f.reader.readChunk(first.fingerprint) == f.frames.last) }
            // Cancel only after one real link (and, optionally, its payload)
            // has been consumed. This uses the production pull walker without
            // relying on sleeps, filesystem notifications, or test hooks.
            withUnsafeCurrentTask { $0?.cancel() }
            return try await traversal.next()
        }
        await expectCancellation { _ = try await task.value }
        #expect(try await f.store.reference("archive-prefix") == reference)
    }

    /// This campaign creates 100,001 real payloads and links plus all immutable
    /// descriptors. Direct fixture writes avoid hundreds of thousands of fsync
    /// calls; production resume and extension still use VivoArtifactStore.
    /// Run explicitly with NUMIVIVO_LONG_ARCHIVE_TESTS=1. It is a functional
    /// regression, not a peak-memory measurement or a crash-durability test.
    /// NUMIVIVO_TEST_ARTIFACTS retains a unique child root and a success receipt
    /// for subsequent verification in processes that did not generate fixtures.
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_LONG_ARCHIVE_TESTS"] == "1"))
    func archiveBeyondOneHundredThousandChunksResumesAndPreservesItsPrefix() async throws {
        let artifactsDirectory = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"]
        let retainArtifacts = artifactsDirectory != nil
        let parent: URL
        if let artifactsDirectory {
            guard !artifactsDirectory.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
                throw VivoArtifactValidationError.invalid("NUMIVIVO_TEST_ARTIFACTS must name a nonempty directory")
            }
            parent = URL(fileURLWithPath: artifactsDirectory, isDirectory: true).standardizedFileURL
            try FileManager.default.createDirectory(at: parent, withIntermediateDirectories: true)
        } else {
            parent = FileManager.default.temporaryDirectory
        }
        let root = parent.appendingPathComponent("numivivo-long-archive-\(UUID().uuidString)")
        // mkdir must claim a previously absent child. Never clean up or reuse
        // a path whose ownership this invocation did not establish exclusively.
        guard Darwin.mkdir(root.path, 0o700) == 0 else {
            throw NSError(domain: NSPOSIXErrorDomain, code: Int(errno), userInfo: [NSFilePathErrorKey: root.path])
        }
        defer { if !retainArtifacts { try? FileManager.default.removeItem(at: root) } }
        print("NUMIVIVO_ARCHIVE_ROOT=\(root.path)")
        let store = try VivoArtifactStore(rootURL: root, createIfNeeded: false)
        var sink = DirectFixtureSink(root: root)
        let system = try VivoCanonicalJSON.fingerprint(Data("long archive system".utf8))
        let configuration = try VivoCanonicalJSON.fingerprint(Data("long archive configuration".utf8))
        let count: UInt64 = 100_001
        let benchmarkCount: UInt64 = 10_000
        let first = frame(0, velocities: false)
        var tail: VivoFingerprint?
        var benchmarkArtifact: VivoStoredArtifact?
        var benchmarkReference: VivoArtifactReference?
        for ordinal in 0..<count {
            try Task.checkCancellation()
            tail = try autoreleasepool {
                let current = frame(ordinal, velocities: false)
                let bytes = try VivoMDTrajectoryChunkCodec.encode([current], systemFingerprint: system,
                    configurationFingerprint: configuration, particleCount: 1,
                    includeVelocities: false, firstFrameOrdinal: ordinal)
                let payload = try sink.put(bytes, kind: "md-trajectory-chunk", mediaType: VivoMDTrajectoryChunkCodec.mediaType)
                let link = VivoMDTrajectoryChunkLink(schema: VivoMDTrajectoryChunkLink.schemaID,
                    systemFingerprint: system, configurationFingerprint: configuration, particleCount: 1,
                    includeVelocities: false, chunkOrdinal: ordinal, firstFrameOrdinal: ordinal,
                    frameCount: 1, previous: tail, payload: payload.fingerprint, payloadBytes: payload.byteCount,
                    firstStep: current.stepIndex, lastStep: current.stepIndex,
                    firstTimePS: current.timePS, lastTimePS: current.timePS)
                return try sink.put(VivoCanonicalJSON.encode(link), kind: "md-trajectory-link",
                                    mediaType: "application/vnd.numivivo.md-trajectory-link+json").fingerprint
            }
            if ordinal + 1 == benchmarkCount {
                let last = frame(ordinal, velocities: false)
                let manifest = VivoMDTrajectoryManifest(schema: VivoMDTrajectoryManifest.schemaID,
                    systemFingerprint: system, configurationFingerprint: configuration, particleCount: 1,
                    includeVelocities: false, frameCount: benchmarkCount, chunkCount: benchmarkCount, tail: tail,
                    firstStep: first.stepIndex, lastStep: last.stepIndex,
                    firstTimePS: first.timePS, lastTimePS: last.timePS, sealed: false)
                let stored = try await storeJSON(manifest, kind: "md-trajectory-manifest", store: store)
                benchmarkArtifact = stored
                benchmarkReference = try await store.setReference("benchmark-prefix-10000", to: stored)
            }
        }
        let last = frame(count - 1, velocities: false)
        let manifest = VivoMDTrajectoryManifest(schema: VivoMDTrajectoryManifest.schemaID,
            systemFingerprint: system, configurationFingerprint: configuration, particleCount: 1,
            includeVelocities: false, frameCount: count, chunkCount: count, tail: tail,
            firstStep: first.stepIndex, lastStep: last.stepIndex,
            firstTimePS: first.timePS, lastTimePS: last.timePS, sealed: false)
        let stored = try await storeJSON(manifest, kind: "md-trajectory-manifest", store: store)
        let reference = try await store.setReference("original-prefix", to: stored)
        let originalBytes = try await store.data(for: stored.fingerprint)
        let benchmark = try #require(benchmarkArtifact)
        let benchmarkReader = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: benchmark.fingerprint)
        let benchmarkValidation = try await benchmarkReader.validate(scope: .allPayloads)
        try #require(benchmarkValidation.indexedChunks == benchmarkCount && benchmarkValidation.verifiedPayloads == benchmarkCount)
        let reader = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: stored.fingerprint)
        var materializingLimitRejected = false
        do { _ = try await reader.index() }
        catch VivoArtifactValidationError.invalid(_) { materializingLimitRejected = true }
        try #require(materializingLimitRejected)
        let originalValidation = try await reader.validate(scope: .allPayloads)
        try #require(originalValidation.indexedChunks == count && originalValidation.verifiedPayloads == count)
        let frameBytes = try VivoMDTrajectoryChunkCodec.frameBytes(particleCount: 1, includeVelocities: false)
        let chunkBytes = UInt64(VivoMDTrajectoryChunkCodec.headerBytes + frameBytes)
        try #require(benchmarkValidation.verifiedPayloadBytes == benchmarkCount * chunkBytes)
        try #require(originalValidation.verifiedPayloadBytes == count * chunkBytes)
        let writer = try await VivoMDTrajectoryArchiveWriter.resume(store: store, manifest: stored.fingerprint)
        let next = frame(count, velocities: false)
        try await writer.append(.init(systemFingerprint: system, configurationFingerprint: configuration,
            stepIndex: next.stepIndex, timePS: next.timePS, positionsNM: next.positionsNM,
            velocitiesNMPerPS: [.zero], periodicCell: nil))
        let extended = try await writer.snapshot()
        let extendedReader = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: extended.fingerprint)
        let extendedTail = try #require(extendedReader.manifest.tail)
        try #require(extendedReader.manifest.frameCount == count + 1 && extendedReader.manifest.chunkCount == count + 1)
        try #require(try await extendedReader.readLink(extendedTail).previous == tail)
        try #require(try await extendedReader.readChunk(extendedTail) == [next])
        let extendedValidation = try await extendedReader.validate(scope: .restart)
        try #require(extendedValidation.indexedChunks == count + 1 && extendedValidation.verifiedPayloads == 1)
        try #require(extendedValidation.verifiedPayloadBytes == chunkBytes)
        let originalManifestPreserved = try await store.data(for: stored.fingerprint) == originalBytes
        let originalReferencePreserved = try await store.reference("original-prefix") == reference
        let expectedBenchmarkReference = try #require(benchmarkReference)
        let benchmarkReferencePreserved = try await store.reference("benchmark-prefix-10000") == expectedBenchmarkReference
        try #require(originalManifestPreserved && originalReferencePreserved && benchmarkReferencePreserved)

        // No receipt exists until every advertised check has passed. Its own
        // publication uses a synced immutable rooted write, unlike bulk fixture
        // generation, and refuses to replace an existing receipt filename.
        try Task.checkCancellation()
        let receipt = LongArchiveValidationReceipt(rootPath: root.path, retained: retainArtifacts,
            completedAt: Date(), benchmark10000: .init(reader: benchmarkReader, validation: benchmarkValidation),
            original100001: .init(reader: reader, validation: originalValidation),
            extended100002: .init(reader: extendedReader, validation: extendedValidation),
            materializingLimitRejected: materializingLimitRejected, originalManifestPreserved: originalManifestPreserved,
            originalReferencePreserved: originalReferencePreserved, benchmarkReferencePreserved: benchmarkReferencePreserved)
        let receiptBytes = try VivoCanonicalJSON.encode(receipt)
        let receiptName = "archive-validation-receipt.json"
        let files = try VivoRootedFileStore(rootURL: root, createIfNeeded: false)
        try #require(try files.writeFile(receiptBytes, relative: receiptName, immutable: true))
        print("NUMIVIVO_ARCHIVE_RECEIPT=\(root.appendingPathComponent(receiptName).path)")
        print(String(decoding: receiptBytes, as: UTF8.self))
    }
}

private struct LongArchiveValidationReceipt: Encodable {
    struct Prefix: Encodable {
        let manifestFingerprint: String
        let frameCount: UInt64
        let chunkCount: UInt64
        let validationScope: String
        let indexedChunks: UInt64
        let verifiedPayloads: UInt64
        let verifiedPayloadBytes: UInt64

        init(reader: VivoMDTrajectoryArchiveReader, validation: VivoMDTrajectoryValidation) {
            manifestFingerprint = reader.manifestFingerprint.hex
            frameCount = reader.manifest.frameCount; chunkCount = reader.manifest.chunkCount
            validationScope = validation.scope.rawValue; indexedChunks = validation.indexedChunks
            verifiedPayloads = validation.verifiedPayloads; verifiedPayloadBytes = validation.verifiedPayloadBytes
        }
    }

    let schema = "numivivo.org/test-evidence/md-trajectory-archive/v1"
    let status = "success"
    let rootPath: String
    let retained: Bool
    let completedAt: Date
    let benchmark10000: Prefix
    let original100001: Prefix
    let extended100002: Prefix
    let materializingLimitRejected: Bool
    let originalManifestPreserved: Bool
    let originalReferencePreserved: Bool
    let benchmarkReferencePreserved: Bool
}

/// A complete content-addressed fixture written before any reader can see it.
/// This bypasses production publication durability, never object hashing or the
/// on-disk format, and is confined to the opt-in test-owned root.
private struct DirectFixtureSink {
    let root: URL
    private var directories: Set<String> = []

    init(root: URL) { self.root = root }

    mutating func put(_ data: Data, kind: String, mediaType: String) throws -> VivoStoredArtifact {
        let fingerprint = try VivoCanonicalJSON.fingerprint(data)
        let hex = fingerprint.hex
        let shard = "\(hex.prefix(2))/\(hex.dropFirst(2).prefix(2))"
        let objectPath = "objects/sha256/\(shard)/\(hex)"
        let record = VivoStoredArtifact(fingerprint: fingerprint, kind: kind, mediaType: mediaType,
            byteCount: UInt64(data.count), objectPath: objectPath, createdAt: Date(timeIntervalSince1970: 0), attributes: [:])
        try write(data, to: objectPath)
        try write(VivoCanonicalJSON.encode(record), to: "descriptors/sha256/\(shard)/\(hex).json")
        return record
    }

    private mutating func write(_ data: Data, to relative: String) throws {
        let path = root.appendingPathComponent(relative)
        let directory = path.deletingLastPathComponent()
        if directories.insert(directory.path).inserted {
            try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true)
        }
        try data.write(to: path, options: .withoutOverwriting)
    }
}
