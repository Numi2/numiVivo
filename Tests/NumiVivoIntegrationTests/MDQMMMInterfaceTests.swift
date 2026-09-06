import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct MDQMMMInterfaceTests {
    private func record<T: Encodable>(_ value: T, _ name: String) throws {
        guard let path = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { return }
        let root = URL(fileURLWithPath: path)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        try VivoCanonicalJSON.encode(value).write(to: root.appendingPathComponent(name + ".json"), options: .atomic)
    }
    private func fixture() throws -> (VivoMolecularStructureDocument, VivoClassicalSystem, VivoMDStateSnapshot, VivoQMMMRegionRequest) {
        let c = try #require(VivoElement.from(symbol: "C")), o = try #require(VivoElement.from(symbol: "O"))
        let h = try #require(VivoElement.from(symbol: "H"))
        let elements = [c,c,o,h,h,o,h,h]
        let atoms = elements.enumerated().map { index, element in
            VivoMolecularAtom(index: UInt32(index), name: "A\(index)", element: element,
                              residueIndex: index < 2 ? 0 : (index < 5 ? 1 : 2))
        }
        let atomPositions: [VivoVector3D] = [.init(0.05,0.50,0.50), .init(1.88,0.48,0.54),
            .init(1.94,0.85,0.50), .init(1.99,0.91,0.50), .init(0.01,0.82,0.50),
            .init(1.00,1.50,1.30), .init(1.06,1.56,1.30), .init(1.07,1.47,1.30)]
        let bonds: [VivoMolecularBond] = [.init(atomA: 0, atomB: 1), .init(atomA: 2, atomB: 3),
            .init(atomA: 2, atomB: 4), .init(atomA: 5, atomB: 6), .init(atomA: 5, atomB: 7)]
        let residues = [VivoMolecularResidue(index: 0, name: "LIG", atomIndices: [0,1]),
                        .init(index: 1, name: "HOH", atomIndices: [2,3,4]), .init(index: 2, name: "HOH", atomIndices: [5,6,7])]
        // Deliberately stale structure cell; the accepted NPT sample uses 2 nm.
        let document = try VivoMolecularStructureDocument(structure: .init(identifier: "synthetic-md-qmmm",
            atoms: atoms, bonds: bonds, residues: residues, conformers: [.init(positionsNM: atomPositions)],
            periodicCell: .init(a: .init(3,0,0), b: .init(0,3,0), c: .init(0,0,3))))
        let order = [4,0,7,2,1,6,3,5]
        let inverse = Dictionary(uniqueKeysWithValues: order.enumerated().map { ($0.element, UInt32($0.offset)) })
        var particles = order.enumerated().map { slot, atom in
            VivoClassicalParticle(index: UInt32(slot), atomIndex: UInt32(atom), typeIdentifier: "physical",
                massDa: elements[atom].atomicNumber == 1 ? 1.008 : 12.0,
                chargeE: elements[atom].atomicNumber == 1 ? 0.5 : 0, sigmaNM: 0.12, epsilonKJPerMol: 0.1)
        }
        particles += [VivoClassicalParticle(index: 8, atomIndex: nil, typeIdentifier: "site", role: .virtualSite,
                                           massDa: 0, chargeE: -1, sigmaNM: 0.08, epsilonKJPerMol: 0.02),
                      .init(index: 9, atomIndex: nil, typeIdentifier: "site", role: .virtualSite,
                            massDa: 0, chargeE: -1, sigmaNM: 0.08, epsilonKJPerMol: 0.02)]
        let sites = [VivoLinearVirtualSite(siteParticle: 8, parentParticles: [inverse[2]!,inverse[3]!,inverse[4]!], weights: [-0.2,0.6,0.6]),
                     .init(siteParticle: 9, parentParticles: [inverse[5]!,inverse[6]!,inverse[7]!], weights: [-0.2,0.6,0.6])]
        let system = VivoClassicalSystem(identifier: "synthetic-md-qmmm", structureFingerprint: document.structureFingerprint,
            particles: particles, linearVirtualSites: sites,
            nonbondedExceptions: [.init(a: inverse[0]!, b: inverse[1]!, coulombScale: 0, lennardJonesScale: 0)])
        var positions = order.map { atomPositions[$0] }
        // Reconstruct the boundary-spanning water from whole parents, then wrap site.
        positions.append(.init(0.012,0.868,0.50))
        positions.append(atomPositions[5] * -0.2 + atomPositions[6] * 0.6 + atomPositions[7] * 0.6)
        let fingerprint = try system.fingerprint()
        let snapshot = VivoMDStateSnapshot(systemFingerprint: fingerprint, configurationFingerprint: fingerprint,
            stepIndex: 91, timePS: 0.182, positionsNM: positions,
            velocitiesNMPerPS: Array(repeating: .zero, count: positions.count),
            periodicCell: .init(a: .init(2,0,0), b: .init(0,2,0), c: .init(0,0,2)))
        return (document, system, snapshot, .init(qmAtomIndices: [0], alphaElectrons: 4, betaElectrons: 3))
    }

    /// Translation/rotation-invariant synthetic center potential. This tests the
    /// interface Jacobians, not the scientific accuracy of an electronic solver.
    private func centerModel(_ system: VivoElectronicSystem) throws -> VivoQMMMElectronicForceResult {
        let centers = system.nuclei.map(\.positionBohr) + system.pointCharges.map(\.positionBohr)
        var energy = 0.0, forces = [SIMD3<Double>](repeating: .zero, count: centers.count)
        for i in centers.indices { for j in 0..<i {
            let d = centers[i] - centers[j], k = 0.001 * Double(i + j + 1)
            energy += 0.5 * k * (d.x*d.x + d.y*d.y + d.z*d.z)
            forces[i] -= d * k; forces[j] += d * k
        } }
        let vectors = forces.map { VivoVector3D($0.x,$0.y,$0.z) }
        return try .init(system: system, energyHartree: energy,
            nucleusForcesHartreePerBohr: Array(vectors.prefix(system.nuclei.count)),
            pointChargeForcesHartreePerBohr: Array(vectors.dropFirst(system.nuclei.count)),
            derivativeMethod: "synthetic-center-pair-harmonic-fixture")
    }

    @Test func periodicSolventPreparationUsesCurrentCellAndPreservesEveryMapping() throws {
        let (document, system, snapshot, request) = try fixture()
        let before = snapshot
        let frame = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: snapshot,
                                                   request: request, solventPolicy: .init(radiusNM: 0.40))
        try frame.validate(document: document, system: system)
        try frame.validateSource(document: document, system: system, snapshot: snapshot)
        try record(frame, "md-qmmm-prepared-frame")
        #expect(snapshot == before)
        #expect(frame.sourceStep == 91 && frame.sourceTimePS == 0.182)
        #expect(frame.cluster.sourcePeriodicCell == snapshot.periodicCell)
        #expect(frame.cluster.atomToParticle == [1,4,3,6,0,7,5,2])
        #expect(frame.cluster.particleImages.map(\.particleIndex) == Array(0..<UInt32(system.particles.count)))
        #expect(frame.request.qmAtomIndices == [0,2,3,4])
        #expect(frame.request.alphaElectrons == 9 && frame.request.betaElectrons == 8)
        #expect(frame.solventPromotion?.promotedResidues.map(\.residueIndex) == [1])
        #expect(frame.solventPromotion?.requiresUnwrappedFiniteClusterForQMMM == false)
        #expect(frame.solventPromotion?.periodicCellUsed == true)
        #expect(frame.region.excludedQMParticles.contains(8))
        #expect(frame.region.electronicSystem.pointCharges.contains { $0.classicalParticleIndex == 9 })
        #expect(!frame.region.electronicSystem.pointCharges.contains { $0.classicalParticleIndex == 8 })
        #expect(frame.cluster.particlePositionsNM.count == snapshot.positionsNM.count)
        #expect(frame.cluster.particleImages[0].latticeImage == [0,0,0])
        #expect(frame.cluster.particleImages[3].latticeImage == [-1,0,0])
        #expect(frame.cluster.particleImages[8].latticeImage == nil)
        #expect(abs(frame.cluster.particlePositionsNM[4].x + 0.12) < 1e-12)
        #expect(abs(frame.cluster.particlePositionsNM[8].x - 0.012) < 1e-12)
        for site in system.linearVirtualSites ?? [] {
            let expected = zip(site.parentParticles, site.weights).reduce(VivoVector3D.zero) {
                $0 + frame.cluster.particlePositionsNM[Int($1.0)] * $1.1
            }
            #expect((expected - frame.cluster.particlePositionsNM[Int(site.siteParticle)]).norm < 1e-12)
        }
        let encoded = try VivoCanonicalJSON.encode(frame)
        let decoded = try VivoCanonicalJSON.decode(VivoMDQMMMPreparedFrame.self, from: encoded)
        #expect(decoded == frame)
        try decoded.validate(document: document, system: system)
    }

    @Test func forceProjectionLJAndVirtualSitesMatchTotalEnergyDifferences() throws {
        let (document, system, snapshot, request) = try fixture()
        let frame = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: snapshot,
                                                   request: request, solventPolicy: .init(radiusNM: 0.40))
        let electronic = try centerModel(frame.region.electronicSystem)
        let result = try VivoQMMMForceMapper.assemble(document: document, system: system, frame: frame, electronic: electronic)
        #expect(result.particleForcesHartreePerBohr[8] == .zero && result.particleForcesHartreePerBohr[9] == .zero)
        #expect(result.netForceHartreePerBohr.norm < 1e-10)
        #expect(result.netTorqueHartree.norm < 1e-9)
        let report = try VivoQMMMForceChecks.check(document: document, system: system, frame: frame, electronic: electronic,
            configuration: .init(absoluteForceTolerance: 1e-7, relativeForceTolerance: 1e-6),
            electronicEnergy: { try centerModel($0).energyHartree })
        try record(report, "md-qmmm-force-check")
        #expect(report.passed && report.completePhysicalAtomCheck)
        #expect(report.components.count == 24 && report.energyEvaluations == 97)
        #expect(report.maximumAbsoluteError < 1e-7)
        let factor = VivoAtomicUnits.hartreeInKJPerMol / VivoAtomicUnits.bohrInNM
        #expect((result.particleForcesKJPerMolPerNM[1] - result.particleForcesHartreePerBohr[1] * factor).norm < 1e-12)
        let lj = try VivoQMMMCompiler.lennardJones(document: document, system: system,
            particlePositionsNM: frame.cluster.particlePositionsNM, request: frame.request)
        let raw = try #require(lj.particleForcesHartreePerBohr)
        #expect(raw.reduce(VivoVector3D.zero, +).norm < 1e-12)
        #expect(raw[8].norm > 0) // QM-owned virtual-site LJ is not discarded.
        let mappedAtomForce = result.atomForcesHartreePerBohr[1]
        #expect(mappedAtomForce == result.particleForcesHartreePerBohr[4])
    }

    @Test func fixedDistanceLinkJacobianConservesForceTorqueAndRadialDerivative() throws {
        let q = VivoVector3D(1,-2,0.5), m = VivoVector3D(3,-1,-0.5)
        let u = (m-q)/(m-q).norm, f = VivoVector3D(0.3,0.4,-0.2), d = 1.1
        let split = try VivoQMMMLinkForceProjection.project(qmPositionBohr: q, mmPositionBohr: m,
                                                           linkDistanceBohr: d, linkForceHartreePerBohr: f)
        #expect((split.qmForceHartreePerBohr + split.mmForceHartreePerBohr - f).norm < 1e-14)
        let torque = q.cross(split.qmForceHartreePerBohr) + m.cross(split.mmForceHartreePerBohr)
        #expect((torque - (q + u*d).cross(f)).norm < 1e-14)
        let radial = try VivoQMMMLinkForceProjection.project(qmPositionBohr: q, mmPositionBohr: m,
                                                            linkDistanceBohr: d, linkForceHartreePerBohr: u*2)
        #expect(radial.mmForceHartreePerBohr.norm < 1e-14)
        #expect((radial.qmForceHartreePerBohr - u*2).norm < 1e-14)
    }

    @Test func forceChecksRejectWrongForcesMissingMMReactionsAndStaleIdentity() throws {
        let (document, system, snapshot, request) = try fixture()
        let frame = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: snapshot, request: request)
        let electronic = try centerModel(frame.region.electronicSystem)
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMElectronicForceResult(system: frame.region.electronicSystem, energyHartree: electronic.energyHartree,
                nucleusForcesHartreePerBohr: electronic.nucleusForcesHartreePerBohr,
                pointChargeForcesHartreePerBohr: [], derivativeMethod: "missing-MM-reactions")
        }
        var moved = frame.region.electronicSystem
        moved.nuclei[0].positionBohr.x += 0.01
        #expect(throws: (any Error).self) { try electronic.validate(system: moved) }
        var wrong = electronic.nucleusForcesHartreePerBohr
        wrong[0] = wrong[0] + .init(0.02,0,0)
        let corrupted = try VivoQMMMElectronicForceResult(system: frame.region.electronicSystem, energyHartree: electronic.energyHartree,
            nucleusForcesHartreePerBohr: wrong, pointChargeForcesHartreePerBohr: electronic.pointChargeForcesHartreePerBohr,
            derivativeMethod: "deliberately-wrong-fixture")
        let failed = try VivoQMMMForceChecks.check(document: document, system: system, frame: frame, electronic: corrupted,
            configuration: .init(atomIndices: [0,1]), electronicEnergy: { try centerModel($0).energyHartree })
        #expect(!failed.passed && !failed.completePhysicalAtomCheck)
        #expect(failed.maximumAbsoluteError > 0.01)
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMForceChecks.check(document: document, system: system, frame: frame, electronic: electronic,
                configuration: .init(maximumEnergyEvaluations: 2), electronicEnergy: { try centerModel($0).energyHartree })
        }
    }

    @Test func energyOnlyBackendIncludesEmbeddingCenterDerivatives() throws {
        let (document, system, snapshot, request) = try fixture()
        let frame = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: snapshot, request: request)
        let numerical = try VivoQMMMElectronicFiniteDifferences.evaluate(system: frame.region.electronicSystem,
            energy: { try centerModel($0).energyHartree })
        let reference = try centerModel(frame.region.electronicSystem)
        #expect(numerical.pointChargeForcesHartreePerBohr.count == reference.pointChargeForcesHartreePerBohr.count)
        for (a,b) in zip(numerical.pointChargeForcesHartreePerBohr, reference.pointChargeForcesHartreePerBohr) {
            #expect((a-b).norm < 1e-7)
        }
    }

    @Test func rejectsStaleSitesPartialSolventAndInsufficientCapacity() throws {
        let (document, system, snapshot, request) = try fixture()
        var stale = snapshot; stale.positionsNM[8].x += 0.03
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: stale, request: request)
        }
        var partial = request; partial.qmAtomIndices = [0,2]
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: snapshot,
                request: partial, solventPolicy: .init(radiusNM: 0.4))
        }
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: snapshot,
                request: request, solventPolicy: .init(radiusNM: 0.4, maximumPromotedAtoms: 2))
        }
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: snapshot,
                request: request, configuration: .init(maximumImagePairs: 1))
        }
    }

    @Test func noncontractiblePeriodicCycleFailsInsteadOfSplittingAMolecule() throws {
        let c = try #require(VivoElement.from(symbol: "C"))
        let atoms = (0..<3).map { VivoMolecularAtom(index: UInt32($0), name: "C\($0)", element: c) }
        let positions: [VivoVector3D] = [.init(0,0,0),.init(0.7,0,0),.init(1.4,0,0)]
        let cell = VivoPeriodicCell(a: .init(2,0,0), b: .init(0,2,0), c: .init(0,0,2))
        let document = try VivoMolecularStructureDocument(structure: .init(identifier: "winding", atoms: atoms,
            bonds: [.init(atomA: 0, atomB: 1),.init(atomA: 1, atomB: 2),.init(atomA: 2, atomB: 0)],
            conformers: [.init(positionsNM: positions)], periodicCell: cell))
        let particles = atoms.map { VivoClassicalParticle(index: $0.index, atomIndex: $0.index, typeIdentifier: "C",
            massDa: 12, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0) }
        let system = VivoClassicalSystem(identifier: "winding", structureFingerprint: document.structureFingerprint, particles: particles)
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMClusterBuilder.reconstruct(document: document, system: system, particlePositionsNM: positions,
                periodicCell: cell, sourceFrameFingerprint: system.fingerprint(), qmAtomIndices: [0])
        }
    }

    @Test func skewCellUsesClosestLatticeImageRatherThanFractionalRounding() throws {
        let c = try #require(VivoElement.from(symbol: "C"))
        let cell = VivoPeriodicCell(a: .init(2,0,0), b: .init(1.8,0.5,0), c: .init(0,0,2))
        let positions: [VivoVector3D] = [.zero, .init(1.862,0.245,0)]
        let atoms = (0..<2).map { VivoMolecularAtom(index: UInt32($0), name: "C\($0)", element: c) }
        let document = try VivoMolecularStructureDocument(structure: .init(identifier: "skew-cell", atoms: atoms,
            bonds: [.init(atomA: 0, atomB: 1)], conformers: [.init(positionsNM: positions)], periodicCell: cell))
        let particles = atoms.map { VivoClassicalParticle(index: $0.index, atomIndex: $0.index, typeIdentifier: "C",
            massDa: 12, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0) }
        let system = VivoClassicalSystem(identifier: "skew-cell", structureFingerprint: document.structureFingerprint, particles: particles)
        let cluster = try VivoQMMMClusterBuilder.reconstruct(document: document, system: system,
            particlePositionsNM: positions, periodicCell: cell, sourceFrameFingerprint: system.fingerprint(), qmAtomIndices: [0,1])
        let found = cluster.particlePositionsNM[1] - cluster.particlePositionsNM[0]
        var best = Double.infinity
        for i in -3...3 { for j in -3...3 { for k in -3...3 {
            best = min(best, (positions[1] - cell.a * Double(i) - cell.b * Double(j) - cell.c * Double(k)).norm)
        } } }
        #expect(abs(found.norm - best) < 1e-12)
        #expect(found.norm < 0.27 && positions[1].norm > 1.8)
        #expect(cluster.particleImages[1].latticeImage == [0,-1,0])
    }

    @Test func nilSampleCellIsNotReplacedByTheImportedCell() throws {
        let (document, system, snapshot, request) = try fixture()
        let periodic = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: snapshot,
            request: request, solventPolicy: .init(radiusNM: 0.4))
        var finite = snapshot
        finite.periodicCell = nil; finite.positionsNM = periodic.cluster.particlePositionsNM
        let frame = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: finite,
            request: request, solventPolicy: .init(radiusNM: 0.4, periodicMinimumImage: false))
        #expect(frame.cluster.sourcePeriodicCell == nil && frame.solventPromotion?.periodicCellUsed == false)
        #expect(frame.request.qmAtomIndices == periodic.request.qmAtomIndices)
        for image in frame.cluster.particleImages where image.structureAtomIndex != nil {
            #expect(image.latticeImage == [0,0,0])
        }
        try frame.validateSource(document: document, system: system, snapshot: finite)
    }

    @Test func molecularPromotionClosesAcrossResidueRecords() throws {
        let c = try #require(VivoElement.from(symbol: "C")), h = try #require(VivoElement.from(symbol: "H"))
        let atoms = [VivoMolecularAtom(index: 0, name: "C", element: c, residueIndex: 0),
            .init(index: 1, name: "H1", element: h, residueIndex: 1), .init(index: 2, name: "H2", element: h, residueIndex: 2)]
        let positions: [VivoVector3D] = [.zero,.init(0.10,0,0),.init(0.174,0,0)]
        let document = try VivoMolecularStructureDocument(structure: .init(identifier: "split-solvent-records", atoms: atoms,
            bonds: [.init(atomA: 1, atomB: 2)], residues: [.init(index: 0, name: "LIG", atomIndices: [0]),
                .init(index: 1, name: "H2", atomIndices: [1]), .init(index: 2, name: "H2", atomIndices: [2])],
            conformers: [.init(positionsNM: positions)]))
        let particles = atoms.map { VivoClassicalParticle(index: $0.index, atomIndex: $0.index, typeIdentifier: $0.element.symbol,
            massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0) }
        let system = VivoClassicalSystem(identifier: "split-solvent-records", structureFingerprint: document.structureFingerprint, particles: particles)
        let request = VivoQMMMRegionRequest(qmAtomIndices: [0], alphaElectrons: 3, betaElectrons: 3)
        let policy = VivoQMMMSolventPromotionPolicy(radiusNM: 0.15, solventResidueNames: ["H2"], distanceAtoms: .allAtoms)
        let result = try VivoQMMMSolventPromotion.promote(document: document, system: system,
            particlePositionsNM: positions, base: request, policy: policy)
        #expect(result.qmAtomIndices == [0,1,2] && result.promotedResidues.count == 2)
        #expect(result.addedAlphaElectrons == 1 && result.addedBetaElectrons == 1)
        var tooSmall = policy; tooSmall.maximumPromotedResidues = 1
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMSolventPromotion.promote(document: document, system: system,
                particlePositionsNM: positions, base: request, policy: tooSmall)
        }
    }

    @Test func acceptedCheckpointAndSourceMappingSurviveRoundTrip() throws {
        let (document, system, snapshot, request) = try fixture()
        let positions = snapshot.positionsNM.map { VivoVector3D(Double(Float($0.x)),Double(Float($0.y)),Double(Float($0.z))) }
        let checkpoint = VivoMDCheckpoint(systemFingerprint: snapshot.systemFingerprint,
            configurationFingerprint: snapshot.configurationFingerprint, acceptedStep: snapshot.stepIndex,
            timePS: snapshot.timePS, positionsNM: positions, velocitiesNMPerPS: snapshot.velocitiesNMPerPS,
            periodicCell: snapshot.periodicCell)
        let frame = try VivoQMMMCompiler.prepareMD(document: document, system: system, checkpoint: checkpoint,
            request: request, solventPolicy: .init(radiusNM: 0.4))
        #expect(frame.sourceCheckpointFingerprint == (try checkpoint.fingerprint()))
        let electronic = try centerModel(frame.region.electronicSystem)
        let forces = try VivoQMMMForceMapper.assemble(document: document, system: system, frame: frame, electronic: electronic)
        #expect(forces.sourceCheckpointFingerprint == frame.sourceCheckpointFingerprint)
        #expect(forces.electronicDerivativeMethod == electronic.derivativeMethod)
        // Tampering with an image can preserve the finite coordinates but cannot
        // pass an audit against the original sampled coordinates.
        var rounded = snapshot; rounded.positionsNM = positions
        var object = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(frame)) as? [String: Any])
        var cluster = try #require(object["cluster"] as? [String: Any])
        var images = try #require(cluster["particleImages"] as? [[String: Any]])
        images[0]["latticeImage"] = [1,0,0]; cluster["particleImages"] = images; object["cluster"] = cluster
        let tampered = try VivoCanonicalJSON.decode(VivoMDQMMMPreparedFrame.self, from: JSONSerialization.data(withJSONObject: object))
        #expect(throws: (any Error).self) { try tampered.validateSource(document: document, system: system, snapshot: rounded) }
    }

    @Test func multipleLinkForcesAccumulateOnTheSharedQMEndpoint() throws {
        let c = try #require(VivoElement.from(symbol: "C"))
        let positions: [VivoVector3D] = [.init(0.1,0,0),.init(1.9,0,0),.init(0.1,1.8,0),.init(1,1,0.5)]
        let atoms = (0..<4).map { VivoMolecularAtom(index: UInt32($0), name: "C\($0)", element: c) }
        let cell = VivoPeriodicCell(a: .init(2,0,0), b: .init(0,2,0), c: .init(0,0,2))
        let document = try VivoMolecularStructureDocument(structure: .init(identifier: "two-links", atoms: atoms,
            bonds: [.init(atomA: 0, atomB: 1),.init(atomA: 0, atomB: 2)],
            conformers: [.init(positionsNM: positions)], periodicCell: cell))
        let particles = atoms.map { VivoClassicalParticle(index: $0.index, atomIndex: $0.index, typeIdentifier: "C",
            massDa: 12, chargeE: $0.index == 3 ? 0.25 : 0, sigmaNM: 0, epsilonKJPerMol: 0) }
        let system = VivoClassicalSystem(identifier: "two-links", structureFingerprint: document.structureFingerprint, particles: particles)
        let fingerprint = try system.fingerprint()
        let snapshot = VivoMDStateSnapshot(systemFingerprint: fingerprint, configurationFingerprint: fingerprint,
            stepIndex: 0, timePS: 0, positionsNM: positions, velocitiesNMPerPS: Array(repeating: .zero, count: 4), periodicCell: cell)
        let frame = try VivoQMMMCompiler.prepareMD(document: document, system: system, snapshot: snapshot,
            request: .init(qmAtomIndices: [0], alphaElectrons: 4, betaElectrons: 4))
        #expect(frame.region.links.count == 2 && frame.region.z1ZeroedMMParticles == [1,2])
        let electronic = try centerModel(frame.region.electronicSystem)
        let report = try VivoQMMMForceChecks.check(document: document, system: system, frame: frame, electronic: electronic,
            electronicEnergy: { try centerModel($0).energyHartree })
        #expect(report.passed && report.completePhysicalAtomCheck)
    }
}
