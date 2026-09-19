import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderStructureSourceTests: XCTestCase {
    private var plan: VivoMolecularInterface.Plan {
        .init(binderChains: ["B"], targetChains: ["T"], conformerID: "model-1")
    }
    private func pdb() -> String {
        var lines: [String] = []
        for (residue, chain) in ["B", "T"].enumerated() {
            for (atom, name) in ["N", "CA", "C", "O"].enumerated() {
                let element = name == "N" ? "N" : name == "O" ? "O" : "C"
                let prefix = "ATOM  " + String(format: "%5d", residue * 4 + atom + 1) + " "
                    + name.padding(toLength: 4, withPad: " ", startingAt: 0) + " GLY " + chain + "   1    "
                lines.append(prefix + String(format: "%8.3f%8.3f%8.3f%6.2f%6.2f",
                    Double(atom), Double(residue) * 4, 0.0, 1.0, 0.0) + "           " + element + "  ")
            }
        }
        return lines.joined(separator: "\n") + "\nEND\n"
    }
    private func cif() -> String {
        var text = """
        data_synthetic
        loop_
        _atom_site.group_PDB
        _atom_site.id
        _atom_site.type_symbol
        _atom_site.label_atom_id
        _atom_site.label_comp_id
        _atom_site.label_asym_id
        _atom_site.label_seq_id
        _atom_site.Cartn_x
        _atom_site.Cartn_y
        _atom_site.Cartn_z
        _atom_site.occupancy

        """
        for (residue, chain) in ["B", "T"].enumerated() {
            for (atom, name) in ["N", "CA", "C", "O"].enumerated() {
                let element = name == "N" ? "N" : name == "O" ? "O" : "C"
                text += "ATOM \(residue * 4 + atom + 1) \(element) \(name) GLY \(chain) 1 \(atom) \(residue * 4) 0 1\n"
            }
        }
        return text + "#\n"
    }
    private func input(_ text: String, format: VivoBinderStructureSources.Format = .pdb,
                       digest: String? = nil, expected: [String: String] = ["T": "G"]) throws -> VivoBinderStructureSources.Input {
        .init(sources: [.init(candidateID: "candidate", target: "T", sourceLabel: "synthetic source fixture",
            format: format, contents: text, sha256: try digest ?? VivoCanonicalJSON.fingerprint(Data(text.utf8)).hex,
            targetChainSequences: expected, interfacePlan: plan)])
    }
    private func dataset(sequence: String = "G") -> VivoBinderBenchmark.Dataset {
        .init(sourceSHA256: String(repeating: "a", count: 64), assay: "synthetic", groupingMethod: "synthetic ids",
            records: [.init(id: "candidate", target: "T", leakageGroup: "g", outcome: .binder,
                rawOutcome: "binder", features: ["ipsae_min_boltz2": 0.5], sourceFields: ["sequence": sequence])])
    }

    func testPDBUsesNativeUnitsAndInterfaceOwner() throws {
        let raw = try input(pdb())
        let parsed = try VivoBinderStructureSources.reconstruct(raw)
        let result = try VivoBinderStructuralFeatures.augment(dataset(), input: parsed)
        XCTAssertEqual(result.interfaces["candidate"]?.minimumDistanceNM ?? -1, 0.4, accuracy: 1e-14)
        XCTAssertEqual(result.interfaces["candidate"]?.contactAtomPairs, 14)
        XCTAssertEqual(parsed.observations[0].structure.atoms.count, 8)
        let bytes = try VivoCanonicalJSON.encode(raw)
        XCTAssertEqual(bytes, try VivoCanonicalJSON.encode(VivoCanonicalJSON.decode(
            VivoBinderStructureSources.Input.self, from: bytes)))
    }
    func testMMCIFUsesSameNativeGeometry() throws {
        let parsed = try VivoBinderStructureSources.reconstruct(input(cif(), format: .mmcif))
        let result = try VivoBinderStructuralFeatures.augment(dataset(), input: parsed)
        XCTAssertEqual(result.interfaces["candidate"]?.minimumDistanceNM ?? -1, 0.4, accuracy: 1e-14)
        XCTAssertEqual(result.interfaces["candidate"]?.contactAtomPairs, 14)
    }
    func testSourceHashMismatchRejected() throws {
        XCTAssertThrowsError(try VivoBinderStructureSources.reconstruct(input(pdb(), digest: String(repeating: "0", count: 64))))
    }
    func testWrongDeclaredTargetRejected() throws {
        XCTAssertThrowsError(try VivoBinderStructureSources.reconstruct(input(pdb(), expected: ["T": "A"])))
        XCTAssertThrowsError(try VivoBinderStructureSources.reconstruct(input(pdb(), expected: ["OTHER": "G"])))
    }
    func testBinderStillMustMatchPublishedSequence() throws {
        let parsed = try VivoBinderStructureSources.reconstruct(input(pdb()))
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(dataset(sequence: "A"), input: parsed))
    }
    func testDuplicateSourceCandidateRejected() throws {
        let source = try input(pdb()).sources[0]
        XCTAssertThrowsError(try VivoBinderStructureSources.reconstruct(.init(sources: [source, source])))
    }
    func testCRLFBytesRemainSourceIdentity() throws {
        let original = pdb(), windows = original.replacingOccurrences(of: "\n", with: "\r\n")
        XCTAssertNotEqual(try input(original).sources[0].sha256, try input(windows).sources[0].sha256)
        let parsed = try VivoBinderStructureSources.reconstruct(input(windows))
        XCTAssertEqual(parsed.observations[0].structure.atoms.count, 8)
    }
    func testSourceNULAndEmptyInputRejected() throws {
        XCTAssertThrowsError(try VivoBinderStructureSources.reconstruct(input(pdb() + "\0")))
        XCTAssertThrowsError(try VivoBinderStructureSources.reconstruct(.init(sources: [])))
    }
}
