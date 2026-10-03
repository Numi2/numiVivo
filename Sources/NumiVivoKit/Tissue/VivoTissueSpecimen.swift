import Foundation

/// Shared biological data contract. Coordinates and identity are independent of
/// assay modality; measured, inferred and hypothetical values never share a flag.
public enum VivoEvidenceState: String, Codable, Sendable {
    case measured = "MEASURED"
    case validatedDomain = "SIMULATED — VALIDATED DOMAIN"
    case outOfDistribution = "SIMULATED — OUT OF DISTRIBUTION"
    case inference = "MODEL INFERENCE"
    case hypothesis = "HYPOTHESIS"
    case unavailable = "UNAVAILABLE"
}
public struct VivoTissueEvidence: Codable, Sendable {
    public let state: VivoEvidenceState
    public let sourceID: String
    public let modelID: String?
    public let assumptions: [String]
    public let validationDomain: String?
    public let uncertainty: String
    public let heldOut: Bool
}
public struct VivoTissueSource: Codable, Sendable {
    public let id: String, uri: String, sha256: String, description: String
}
public struct VivoTissueEntity: Codable, Sendable {
    public let id: String, level: String
    public let parentID: String?
    public let memberIDs: [String]
    public let position: [Double]?
    public let biologicalUnitID: String
    public let annotations: [String:String]
    public let evidence: VivoTissueEvidence
}
public struct VivoTissueMeasurement: Codable, Sendable {
    public let id: String, modality: String, unit: String, timepoint: String
    public let entityIDs: [String], featureIDs: [String]
    /// Sparse values: absence is missing unless explicitly declared as zero.
    public let rowOffsets: [Int], featureIndices: [Int], values: [Double]
    public let absentValue: String
    public let evidence: VivoTissueEvidence
}
public struct VivoTissueCoordinateSystem: Codable, Sendable {
    public let id: String, unit: String, axes: [String], description: String
    public let micrometresPerUnit: Double?
}
public struct VivoTissueSpecimen: Codable, Sendable {
    public let format: String, id: String, title: String, organism: String
    public let coordinateSystem: VivoTissueCoordinateSystem
    public let sources: [VivoTissueSource]
    public let entities: [VivoTissueEntity]
    public let measurements: [VivoTissueMeasurement]
    public let limitations: [String]

    public func validate() throws {
        func need(_ b: Bool,_ s:String) throws { if !b { throw NSError(domain:"NumiVivo.Tissue",code:1,userInfo:[NSLocalizedDescriptionKey:s]) } }
        try need(format == "numivivo-tissue-specimen/v1" && !id.isEmpty && !organism.isEmpty,"specimen identity/schema")
        try need(!sources.isEmpty && Set(sources.map(\.id)).count == sources.count,"source identities")
        let sourceIDs=Set(sources.map(\.id))
        for s in sources { try need(s.sha256.count==64 && s.sha256.allSatisfy{$0.isHexDigit} && !s.uri.isEmpty,"source digest/provenance") }
        try need((2...3).contains(coordinateSystem.axes.count) && !coordinateSystem.unit.isEmpty && Set(coordinateSystem.axes).count==coordinateSystem.axes.count,"coordinate axes/unit")
        if let scale=coordinateSystem.micrometresPerUnit { try need(scale.isFinite && scale>0,"coordinate scale") }
        try need(!entities.isEmpty && entities.count<=1_000_000,"entity budget")
        let ids=Set(entities.map(\.id));try need(ids.count==entities.count,"duplicate entity")
        let lookup=Dictionary(uniqueKeysWithValues:entities.map{($0.id,$0)})
        func evidence(_ e:VivoTissueEvidence) throws {
            try need(sourceIDs.contains(e.sourceID) && !e.uncertainty.isEmpty,"evidence provenance/uncertainty")
            if [.validatedDomain,.outOfDistribution,.inference].contains(e.state) { try need(e.modelID != nil,"model evidence needs identity") }
            if e.state == .validatedDomain { try need(e.validationDomain != nil,"validated output needs an explicit domain") }
        }
        let levels=["molecule","cell","neighborhood","region","specimen"]
        for e in entities {
            try need(!e.id.isEmpty && levels.contains(e.level) && !e.biologicalUnitID.isEmpty,"entity identity/level/unit")
            if let p=e.parentID { try need(ids.contains(p) && p != e.id,"entity parent")
                var visited=Set([e.id]);var next:String?=p
                while let n=next { try need(visited.insert(n).inserted,"cyclic hierarchy");next=lookup[n]?.parentID }
            }
            try need(Set(e.memberIDs).count==e.memberIDs.count && e.memberIDs.allSatisfy{ids.contains($0) && $0 != e.id},"neighborhood membership")
            if let xyz=e.position { try need(xyz.count==coordinateSystem.axes.count && xyz.allSatisfy{$0.isFinite},"coordinate dimensions") }
            try evidence(e.evidence)
        }
        try need(Set(measurements.map(\.id)).count==measurements.count,"measurement identity")
        for m in measurements {
            try need(["rna","protein","morphology","extracellular-field","image-reference"].contains(m.modality) && !m.unit.isEmpty && !m.timepoint.isEmpty,"measurement modality/units/time")
            try need(!m.entityIDs.isEmpty && !m.featureIDs.isEmpty && Set(m.entityIDs).count==m.entityIDs.count && Set(m.featureIDs).count==m.featureIDs.count,"measurement axes")
            try need(m.entityIDs.allSatisfy{ids.contains($0)} && ["zero","missing"].contains(m.absentValue),"measurement identities/missingness")
            try need(m.rowOffsets.count==m.entityIDs.count+1 && m.rowOffsets.first==0 && m.rowOffsets.last==m.values.count && m.values.count==m.featureIndices.count && m.values.count<=20_000_000,"sparse shape")
            try need(m.values.allSatisfy{$0.isFinite} && m.featureIndices.allSatisfy{m.featureIDs.indices.contains($0)},"sparse values")
            for i in m.entityIDs.indices {
                let a=m.rowOffsets[i],b=m.rowOffsets[i+1];try need(a>=0 && b>=a && b<=m.values.count,"sparse offsets")
                try need(zip(m.featureIndices[a..<b],m.featureIndices[a..<b].dropFirst()).allSatisfy{$0<$1},"duplicate/unsorted feature index")
            }
            if m.unit=="umiCount" { try need(m.modality=="rna" && m.values.allSatisfy{$0>=0 && $0.rounded()==$0},"RNA integer counts") }
            try evidence(m.evidence)
        }
    }
}
