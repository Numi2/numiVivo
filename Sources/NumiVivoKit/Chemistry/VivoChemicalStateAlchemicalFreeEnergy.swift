import Foundation

public struct VivoChemicalAlchemicalState: Codable, Sendable, Equatable {
    public var identifier: String
    public var boundProtonOffset: Int

    public init(identifier: String, boundProtonOffset: Int) {
        self.identifier = identifier
        self.boundProtonOffset = boundProtonOffset
    }

    public func validate() throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 512,
              (-128...128).contains(boundProtonOffset) else {
            throw VivoChemistryError.invalid("alchemical chemical-state identity")
        }
    }
}

public struct VivoProtonReservoirCalibration: Codable, Sendable, Equatable {
    /// pH at which this calibration is defined.
    public var referencePH: Double
    /// Additive semigrand contribution per relative bound proton at referencePH.
    /// If state i has n_i more bound protons than state zero, n_i*this value is
    /// added to its physical relative free energy before pH reweighting.
    public var contributionPerBoundProtonKJPerMol: Double
    /// Optional uncertainty of the per-proton contribution. Independence from
    /// MBAR sampling uncertainty is an explicit approximation in the reported SD.
    public var standardDeviationKJPerMol: Double?
    public var origin: VivoKineticOrigin
    public var evidence: VivoKineticEvidence

    public init(referencePH: Double,
                contributionPerBoundProtonKJPerMol: Double,
                standardDeviationKJPerMol: Double? = nil,
                origin: VivoKineticOrigin,
                evidence: VivoKineticEvidence) {
        self.referencePH = referencePH
        self.contributionPerBoundProtonKJPerMol = contributionPerBoundProtonKJPerMol
        self.standardDeviationKJPerMol = standardDeviationKJPerMol
        self.origin = origin
        self.evidence = evidence
    }

    public func validate() throws {
        guard referencePH.isFinite, (-10...30).contains(referencePH),
              contributionPerBoundProtonKJPerMol.isFinite else {
            throw VivoChemistryError.invalid("proton-reservoir calibration")
        }
        if let standardDeviationKJPerMol {
            guard standardDeviationKJPerMol.isFinite, standardDeviationKJPerMol >= 0 else {
                throw VivoChemistryError.invalid("proton-reservoir calibration uncertainty")
            }
        }
        try evidence.validate(origin: origin)
    }
}

public struct VivoChemicalStateAlchemicalFreeEnergyRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/chemical-state-alchemical-free-energy/v1"
    public var schema: String
    public var identifier: String
    public var temperatureK: Double
    public var referencePH: Double
    public var states: [VivoChemicalAlchemicalState]
    public var samples: [VivoMBARReducedPotentialSample]
    public var mbar: VivoMBARConfiguration
    public var protonReservoir: VivoProtonReservoirCalibration?
    public var samplingOrigin: VivoKineticOrigin
    public var samplingEvidence: VivoKineticEvidence

    public init(identifier: String, temperatureK: Double, referencePH: Double,
                states: [VivoChemicalAlchemicalState],
                samples: [VivoMBARReducedPotentialSample],
                mbar: VivoMBARConfiguration = .init(),
                protonReservoir: VivoProtonReservoirCalibration? = nil,
                samplingOrigin: VivoKineticOrigin,
                samplingEvidence: VivoKineticEvidence) {
        schema = Self.schema
        self.identifier = identifier
        self.temperatureK = temperatureK
        self.referencePH = referencePH
        self.states = states
        self.samples = samples
        self.mbar = mbar
        self.protonReservoir = protonReservoir
        self.samplingOrigin = samplingOrigin
        self.samplingEvidence = samplingEvidence
    }
}

public struct VivoChemicalStateAlchemicalEstimate: Codable, Sendable, Equatable {
    public let identifier: String
    public let boundProtonOffset: Int
    /// Physical free energy from MBAR before proton-reservoir calibration.
    public let relativePhysicalFreeEnergyKJPerMol: Double
    public let relativeSemigrandFreeEnergyKJPerMolAtReferencePH: Double
    /// Conditional finite-sample + calibration SD. Model/state completeness and
    /// covariance between the two components are excluded.
    public let conditionalStandardDeviationKJPerMol: Double?
    public let effectiveSamples: Double
}

public struct VivoChemicalStateAlchemicalFreeEnergyResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/chemical-state-alchemical-free-energy-result/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let mbar: VivoMBARResult
    public let estimates: [VivoChemicalStateAlchemicalEstimate]
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint

    public func thermodynamicsRequest(identifier: String,
                                      temperatureK: Double,
                                      referencePH: Double,
                                      targetPH: Double) throws -> VivoQMMMChemicalStateThermodynamicsRequest {
        guard temperatureK.isFinite, referencePH.isFinite, targetPH.isFinite else {
            throw VivoChemistryError.invalid("alchemical thermodynamics adapter context")
        }
        let evidence = VivoKineticEvidence(source: "NumiVivo alchemical multistate free-energy calculation",
            locator: "MBAR relative physical free energies plus explicit proton-reservoir calibration",
            sourceFingerprint: evidenceFingerprint.hex)
        let states = estimates.map { estimate in
            VivoQMMMChemicalThermodynamicState(identifier: estimate.identifier,
                boundProtonOffset: estimate.boundProtonOffset,
                relativeSemigrandFreeEnergyKJPerMol: estimate.relativeSemigrandFreeEnergyKJPerMolAtReferencePH,
                origin: .calculated, evidence: evidence)
        }
        return .init(identifier: identifier, temperatureK: temperatureK,
                     referencePH: referencePH, targetPH: targetPH, states: states)
    }
}

/// Converts a general physical-state MBAR result into semigrand state free
/// energies. Proton-changing states require an explicit reservoir calibration;
/// the MBAR result alone is never interpreted as a pKa or proton chemical potential.
public enum VivoChemicalStateAlchemicalFreeEnergy {
    private static let gasConstantKJ = 0.00831446261815324
    public static let interpretation = "Relative physical chemical-state free energies from decorrelated multistate reduced-potential samples, promoted to reference-pH semigrand free energies only through an explicit evidence-backed proton-reservoir calibration when proton stoichiometry differs. Finite-sample bootstrap uncertainty is retained; state completeness, force-field/electronic accuracy, finite-size corrections and proton-reservoir model error beyond the supplied uncertainty remain separate."

    public static func calculate(_ request: VivoChemicalStateAlchemicalFreeEnergyRequest) throws -> VivoChemicalStateAlchemicalFreeEnergyResult {
        guard request.schema == VivoChemicalStateAlchemicalFreeEnergyRequest.schema,
              !request.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              request.identifier.utf8.count <= 512,
              request.temperatureK.isFinite, request.temperatureK > 0,
              request.referencePH.isFinite, (-10...30).contains(request.referencePH),
              request.states.count >= 2, request.states.count <= 4096,
              Set(request.states.map(\.identifier)).count == request.states.count else {
            throw VivoChemistryError.invalid("chemical-state alchemical free-energy request")
        }
        for state in request.states { try state.validate() }
        try request.mbar.validate()
        try request.samplingEvidence.validate(origin: request.samplingOrigin)
        let protonOffsets = Set(request.states.map(\.boundProtonOffset))
        if protonOffsets.count > 1 {
            guard let calibration = request.protonReservoir else {
                throw VivoChemistryError.unsupported("proton-changing alchemical states require an explicit proton-reservoir calibration")
            }
            try calibration.validate()
            guard calibration.referencePH == request.referencePH else {
                throw VivoChemistryError.invalid("proton-reservoir calibration pH differs from alchemical reference pH")
            }
        } else if let calibration = request.protonReservoir {
            try calibration.validate()
            guard calibration.referencePH == request.referencePH else {
                throw VivoChemistryError.invalid("optional proton-reservoir calibration pH differs from alchemical reference pH")
            }
        }
        let mbar = try VivoMultistateMBAR.solve(samples: request.samples,
                                                stateCount: request.states.count,
                                                configuration: request.mbar)
        guard mbar.converged else {
            throw VivoChemistryError.convergence("chemical-state MBAR did not pass residual, overlap or effective-sample gates")
        }
        let rt = Self.gasConstantKJ * request.temperatureK
        let referenceProtons = request.states[0].boundProtonOffset
        let calibrationValue = request.protonReservoir?.contributionPerBoundProtonKJPerMol ?? 0
        let calibrationSD = request.protonReservoir?.standardDeviationKJPerMol
        var rawSemigrand: [Double] = []
        rawSemigrand.reserveCapacity(request.states.count)
        for index in request.states.indices {
            let physical = mbar.estimates[index].relativeFreeEnergy * rt
            let deltaProtons = request.states[index].boundProtonOffset - referenceProtons
            rawSemigrand.append(physical + Double(deltaProtons) * calibrationValue)
        }
        let referenceSemigrand = rawSemigrand[0]
        let estimates = request.states.indices.map { index -> VivoChemicalStateAlchemicalEstimate in
            let physical = mbar.estimates[index].relativeFreeEnergy * rt
            let deltaProtons = request.states[index].boundProtonOffset - referenceProtons
            let sampleSD = mbar.estimates[index].bootstrapStandardDeviation.map { $0 * rt }
            let combinedSD: Double?
            if sampleSD != nil || (calibrationSD != nil && deltaProtons != 0) {
                let a = sampleSD ?? 0
                let b = Double(abs(deltaProtons)) * (calibrationSD ?? 0)
                combinedSD = sqrt(a * a + b * b)
            } else { combinedSD = nil }
            return .init(identifier: request.states[index].identifier,
                boundProtonOffset: request.states[index].boundProtonOffset,
                relativePhysicalFreeEnergyKJPerMol: physical,
                relativeSemigrandFreeEnergyKJPerMolAtReferencePH: rawSemigrand[index] - referenceSemigrand,
                conditionalStandardDeviationKJPerMol: combinedSD,
                effectiveSamples: mbar.estimates[index].effectiveSamples)
        }
        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Codable {
            let schema: String
            let request: VivoChemicalStateAlchemicalFreeEnergyRequest
            let mbar: VivoMBARResult
            let estimates: [VivoChemicalStateAlchemicalEstimate]
        }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema: "numivivo.org/chemical-state-alchemical-free-energy-evidence/v1",
            request: request, mbar: mbar, estimates: estimates)))
        return .init(schema: VivoChemicalStateAlchemicalFreeEnergyResult.schema,
            requestFingerprint: requestID, mbar: mbar, estimates: estimates,
            interpretation: interpretation, evidenceFingerprint: evidenceID)
    }

    public static func validate(_ result: VivoChemicalStateAlchemicalFreeEnergyResult,
                                request: VivoChemicalStateAlchemicalFreeEnergyRequest) throws {
        guard result == (try calculate(request)) else {
            throw VivoChemistryError.invalid("chemical-state alchemical free energy does not reconstruct")
        }
    }
}