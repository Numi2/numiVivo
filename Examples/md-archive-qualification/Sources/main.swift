import Foundation
import NumiVivoKit

struct Failure: Error { let message: String }
var checks = 0
func require(_ condition: @autoclosure () throws -> Bool, _ label: String) throws {
    guard try condition() else { throw Failure(message: label) }
    checks += 1
}
func rejects(_ label: String, _ operation: () throws -> Void) throws {
    var rejected = false
    do { try operation() } catch { rejected = true }
    try require(rejected, label)
}

let directory = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-md-archive-\(UUID().uuidString)")
try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true)
defer { try? FileManager.default.removeItem(at: directory) }
let systemIdentity = try VivoCanonicalJSON.fingerprint(Data("synthetic-system".utf8))
let configurationIdentity = try VivoCanonicalJSON.fingerprint(Data("synthetic-config".utf8))
let positions: [VivoVector3D] = [.init(0, 0, 0), .init(0.125, 0.25, 0.5)]
let velocities: [VivoVector3D] = [.init(0.5, -0.25, 0.125), .init(-0.5, 0.25, -0.125)]
func snapshot(step: UInt64, time: Double, positions: [VivoVector3D]) -> VivoMDStateSnapshot {
    .init(systemFingerprint: systemIdentity, configurationFingerprint: configurationIdentity,
          stepIndex: step, timePS: time, positionsNM: positions,
          velocitiesNMPerPS: velocities, periodicCell: nil)
}

let url = directory.appendingPathComponent("roundtrip.nvmd")
let header = VivoMDTrajectoryArchiveHeader(systemFingerprint: systemIdentity,
                                          particleCount: 2, includesVelocities: true)
let writer = try VivoMDTrajectoryArchiveWriter(url: url, header: header)
try writer.append(snapshot(step: 0, time: 0, positions: positions), stageIdentifier: "nvt")
try writer.append(snapshot(step: 1, time: 0.001, positions: positions), stageIdentifier: "nvt")
let receipt = try writer.finish()
try require(receipt.frameCount == 2, "writer frame count")
let reader = try VivoMDTrajectoryArchiveReader(url: url)
let first = try reader.next()
try require(first?.positionsNM == positions, "lossless positions")
try require(first?.velocitiesNMPerPS == velocities, "lossless velocities")
try require(reader.verifiedReceipt() == nil, "partial reading is not complete verification")
let second = try reader.next()
try require(second?.metadata.stepIndex == 1, "step metadata")
try require(try reader.next() == nil, "verified footer EOF")
try require(reader.verifiedReceipt() == receipt, "whole-file digest and chain receipt")
try rejects("exclusive publication refuses existing destination") {
    _ = try VivoMDTrajectoryArchiveWriter(url: url, header: header)
}

let bytes = try Data(contentsOf: url)
let headerLength = bytes.withUnsafeBytes {
    Int(UInt32(littleEndian: $0.loadUnaligned(fromByteOffset: 8, as: UInt32.self)))
}
var corrupted = bytes
let firstFrame = 12 + headerLength + 32
corrupted[firstFrame + 16] ^= 1
let corruptURL = directory.appendingPathComponent("corrupt.nvmd")
try corrupted.write(to: corruptURL)
try rejects("frame corruption is rejected") {
    let input = try VivoMDTrajectoryArchiveReader(url: corruptURL)
    _ = try input.next()
}
let truncatedURL = directory.appendingPathComponent("truncated.nvmd")
try Data(bytes.dropLast()).write(to: truncatedURL)
try rejects("truncated footer is not successful EOF") {
    let input = try VivoMDTrajectoryArchiveReader(url: truncatedURL)
    while try input.next() != nil {}
}
var extra = bytes; extra.append(0)
let extraURL = directory.appendingPathComponent("trailing.nvmd")
try extra.write(to: extraURL)
try rejects("trailing bytes are rejected") {
    let input = try VivoMDTrajectoryArchiveReader(url: extraURL)
    while try input.next() != nil {}
}

let emptyURL = directory.appendingPathComponent("empty.nvmd")
let emptyWriter = try VivoMDTrajectoryArchiveWriter(url: emptyURL, header: header)
let emptyReceipt = try emptyWriter.finish()
let emptyReader = try VivoMDTrajectoryArchiveReader(url: emptyURL)
try require(try emptyReader.next() == nil, "zero-frame archive footer")
try require(emptyReader.verifiedReceipt() == emptyReceipt, "zero-frame receipt")

let lossyURL = directory.appendingPathComponent("lossy.nvmd")
let lossy = try VivoMDTrajectoryArchiveWriter(url: lossyURL, header: header)
try rejects("arbitrary Double coordinates are not silently quantized") {
    try lossy.append(snapshot(step: 0, time: 0, positions: [.init(0.1, 0, 0), positions[1]]))
}
lossy.abort()
try require(!FileManager.default.fileExists(atPath: lossyURL.path), "failed archive is not published")

let particles = [
    VivoClassicalParticle(index: 0, atomIndex: 0, typeIdentifier: "H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0),
    VivoClassicalParticle(index: 1, atomIndex: 1, typeIdentifier: "H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
]
let system = VivoClassicalSystem(identifier: "thermal-projection-fixture", structureFingerprint: systemIdentity,
                                 particles: particles,
                                 constraints: [.init(a: 0, b: 1, distanceNM: 0.125)])
let recipe = VivoMDVelocityInitialization(temperatureK: 300, seed: 42, removeCenterOfMass: true)
let coordinates: [VivoVector3D] = [.zero, .init(0.125, 0, 0)]
let thermal = try VivoMDVelocityInitializer.make(system: system, positionsNM: coordinates,
                                                periodicCell: nil, recipe: recipe)
let repeated = try VivoMDVelocityInitializer.make(system: system, positionsNM: coordinates,
                                                 periodicCell: nil, recipe: recipe)
try require(thermal.velocitiesNMPerPS == repeated.velocitiesNMPerPS, "seeded preparation reproducibility")
let v0 = thermal.velocitiesNMPerPS[0], v1 = thermal.velocitiesNMPerPS[1]
try require((v0 + v1).norm < 1e-10, "COM momentum projection")
try require(abs(v0.x - v1.x) < 1e-10, "constraint tangent projection")
try require(thermal.maximumNormalizedConstraintResidual <= 1e-8, "reported projection residual")

let definition = try VivoMDProtocol.biomolecular(identifier: "protocol-fixture", nvtSteps: 5, nptSteps: 5, productionSteps: 5)
try require(definition.stages.count == 4, "four-stage protocol")
let seeds = definition.stages.filter { $0.kind == .dynamics }.map { $0.configuration.randomSeed }
try require(Set(seeds).count == seeds.count, "distinct stage RNG namespaces")
try require(definition.stages.filter { $0.initializeVelocities != nil }.count == 1, "one thermal initialization")

print("\(checks) archive, protocol and preparation assertions passed.")
print("No MD integration, Metal force execution, PME or ensemble qualification was performed by this program.")
