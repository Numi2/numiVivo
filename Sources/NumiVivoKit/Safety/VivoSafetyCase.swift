import Foundation

public enum VivoAssuranceStatus: String, Codable, Sendable, CaseIterable {
    case satisfied, conditional, unsupported, notApplicable, notEvaluated
}

public enum VivoSafetyDisposition: String, Codable, Sendable {
    case computationalResearchOnly, conditionalComputationalUse, rejected
}

public struct VivoSafetyClaim: Codable, Sendable, Equatable, Hashable {
    public let identifier: String
    public let title: String
    public let status: VivoAssuranceStatus
    public let rationale: String
    public let evidence: [String]
    public let assumptions: [String]

    public init(
        identifier: String,
        title: String,
        status: VivoAssuranceStatus,
        rationale: String,
        evidence: [String] = [],
        assumptions: [String] = []
    ) {
        self.identifier = identifier
        self.title = title
        self.status = status
        self.rationale = rationale
        self.evidence = evidence
        self.assumptions = assumptions
    }
}

public struct VivoSafetyFinding: Codable, Sendable, Equatable, Hashable {
    public enum Severity: String, Codable, Sendable, CaseIterable {
        case note, warning, error, blocking
    }

    public let identifier: String
    public let severity: Severity
    public let category: String
    public let subject: String
    public let message: String
    public let requiredMitigation: String

    public init(
        identifier: String,
        severity: Severity,
        category: String,
        subject: String,
        message: String,
        requiredMitigation: String
    ) {
        self.identifier = identifier
        self.severity = severity
        self.category = category
        self.subject = subject
        self.message = message
        self.requiredMitigation = requiredMitigation
    }
}

public struct VivoSafetyCoverage: Codable, Sendable, Equatable {
    public let outputSpecies: [String]
    public let monitoredOutputSpecies: [String]
    public let measuredOutputSpecies: [String]
    public let terminationMonitors: [String]
    public let requiredCouplingChannels: [String]
    public let transactionalRequiredCouplingChannels: [String]
    public let parameterEvidenceCounts: [String: UInt32]
}

public struct VivoSafetyCase: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/safety-case/v1"

    public let schema: String
    public let generatedAt: Date
    public let programFingerprint: VivoFingerprint
    public let programSourceFingerprint: VivoFingerprint
    public let hostContextFingerprint: VivoFingerprint?
    public let experimentFingerprint: VivoFingerprint?
    public let couplingFingerprint: VivoFingerprint?
    public let disposition: VivoSafetyDisposition
    public let claims: [VivoSafetyClaim]
    public let findings: [VivoSafetyFinding]
    public let coverage: VivoSafetyCoverage
    public let limitations: [String]

    public init(
        generatedAt: Date = Date(),
        programFingerprint: VivoFingerprint,
        programSourceFingerprint: VivoFingerprint,
        hostContextFingerprint: VivoFingerprint?,
        experimentFingerprint: VivoFingerprint?,
        couplingFingerprint: VivoFingerprint?,
        disposition: VivoSafetyDisposition,
        claims: [VivoSafetyClaim],
        findings: [VivoSafetyFinding],
        coverage: VivoSafetyCoverage,
        limitations: [String]
    ) {
        self.schema = Self.schema
        self.generatedAt = generatedAt
        self.programFingerprint = programFingerprint
        self.programSourceFingerprint = programSourceFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.experimentFingerprint = experimentFingerprint
        self.couplingFingerprint = couplingFingerprint
        self.disposition = disposition
        self.claims = claims
        self.findings = findings
        self.coverage = coverage
        self.limitations = limitations
    }

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }

    public func envelope(
        metadata: VivoArtifactMetadata,
        createdBy: String,
        softwareVersion: String
    ) -> VivoArtifactEnvelope<VivoSafetyCase> {
        VivoArtifactEnvelope(
            kind: .safetyCase,
            metadata: metadata,
            provenance: .init(
                createdBy: createdBy,
                softwareVersion: softwareVersion,
                sourceArtifacts: [programFingerprint] +
                    [hostContextFingerprint, experimentFingerprint, couplingFingerprint].compactMap { $0 }
            ),
            payload: self
        )
    }
}

public struct VivoSafetyCaseBuilder: Sendable {
    public init() {}

    public func build(
        programPack: VivoProgramPack,
        hostContext: VivoHostContextPack? = nil,
        experiment: VivoExperimentPack? = nil,
        coupling: VivoCouplingPack? = nil,
        generatedAt: Date = Date()
    ) throws -> VivoSafetyCase {
        let species = try programPack.speciesMetadata()
        let parameters = try programPack.parameterMetadata()
        let actions = try programPack.actionMetadata()
        let monitors = try programPack.monitorMetadata()
        let instructions = try programPack.expressionInstructions()
        let hostFingerprint = try fingerprint(hostContext)
        let experimentFingerprint = try fingerprint(experiment)
        let couplingFingerprint = try fingerprint(coupling)

        var claims: [VivoSafetyClaim] = [.init(
            identifier: "SC-PACK-001",
            title: "The evaluated program is an immutable parsed ProgramPack",
            status: .satisfied,
            rationale: "VivoProgramPack construction validates binary identity, section bounds, ABI strides, required sections, and runtime counts.",
            evidence: [
                "program:\(programPack.header.contentFingerprint.hex)",
                "source:\(programPack.header.sourceFingerprint.hex)",
                "compilerABI:\(programPack.header.compilerABI)"
            ]
        )]
        var findings: [VivoSafetyFinding] = []
        var limitations = [
            "This safety case evaluates computational contracts and declared evidence only.",
            "It does not establish biological realizability, efficacy, clinical safety, environmental safety, manufacturing quality, or regulatory authorization.",
            "An absence of detected model violations is not evidence that unmodelled failure modes are absent."
        ]

        let outputs = species.enumerated().filter { $0.element.isOutput }
        let outputNames = outputs.map { $0.element.identifier }.sorted()
        let monitoredIndices = monitorReferences(
            monitors: monitors,
            instructions: instructions,
            speciesCount: species.count,
            findings: &findings
        )
        let monitoredOutputs = outputs
            .filter { monitoredIndices.contains(UInt32($0.offset)) }
            .map { $0.element.identifier }
            .sorted()
        let unmonitored = outputNames.filter { !monitoredOutputs.contains($0) }
        claims.append(.init(
            identifier: "SC-MON-001",
            title: "Declared output species are referenced by runtime monitors",
            status: outputNames.isEmpty ? .notApplicable : (unmonitored.isEmpty ? .satisfied : .unsupported),
            rationale: "\(monitoredOutputs.count) of \(outputNames.count) output species have direct monitor coverage.",
            evidence: monitoredOutputs,
            assumptions: unmonitored
        ))
        unmonitored.forEach {
            add(&findings, "SC-MON-MISSING-\($0)", .blocking, "monitor-coverage", $0,
                "Output species is not loaded by any runtime monitor.",
                "Add rate, state, cumulative-exposure, or downstream-effect constraints that reference this output.")
        }

        let termination = monitors.filter { $0.isTerminationRule }.map { $0.identifier }.sorted()
        claims.append(.init(
            identifier: "SC-TERM-001",
            title: "The program has a declared terminal response",
            status: termination.isEmpty ? .unsupported : .satisfied,
            rationale: termination.isEmpty ? "No termination monitor is compiled." : "Compiled termination monitors are present.",
            evidence: termination
        ))
        if termination.isEmpty {
            add(&findings, "SC-TERM-MISSING", .blocking, "termination", "program",
                "ProgramPack contains no termination monitor.",
                "Declare an independently triggered reversible or permanent shutdown condition.")
        }

        let boundedOutputs = outputs.compactMap { pair -> String? in
            let value = pair.element
            let lower = value.minimum.isFinite && value.minimum > -Float.greatestFiniteMagnitude
            let upper = value.maximum.isFinite && value.maximum < Float.greatestFiniteMagnitude
            return lower && upper ? value.identifier : nil
        }.sorted()
        claims.append(.init(
            identifier: "SC-BOUND-001",
            title: "Output state is bounded",
            status: outputNames.isEmpty ? .notApplicable : (boundedOutputs.count == outputNames.count ? .satisfied : .conditional),
            rationale: "\(boundedOutputs.count) of \(outputNames.count) outputs have finite lower and upper bounds.",
            evidence: boundedOutputs
        ))
        for name in outputNames where !boundedOutputs.contains(name) {
            add(&findings, "SC-BOUND-\(name)", .error, "output-bound", name,
                "Output species lacks complete finite bounds.",
                "Declare biologically and numerically justified lower and upper bounds.")
        }

        for action in actions where (action.kind == .addOutput || action.kind == .express) && action.maximumRate <= 0 {
            let subject = action.targetIndex < UInt32(species.count)
                ? species[Int(action.targetIndex)].identifier
                : "action:\(action.index)"
            add(&findings, "SC-RATE-\(action.index)", .error, "unbounded-rate", subject,
                "Accumulating action has no positive maximum-rate contract.",
                "Compile a finite maximumRate and a cumulative output monitor.")
        }

        let irreversibleActions = actions.filter {
            $0.kind == .requestDifferentiation || $0.kind == .permanentShutdown
        }
        let irreversibleProgram = programPack.runtimeContract.featureFlags & (1 << 5) != 0 || !irreversibleActions.isEmpty
        claims.append(.init(
            identifier: "SC-IRREV-001",
            title: "Irreversible program behavior is declared",
            status: irreversibleProgram ? .conditional : .notApplicable,
            rationale: irreversibleProgram
                ? "Irreversible behavior is declared and requires host-context containment evaluation."
                : "No irreversible feature flag or action was found.",
            evidence: irreversibleActions.map { "action:\($0.index):\($0.kind.rawValue)" }
        ))

        let parameterCoverage = evaluateParameters(parameters, claims: &claims, findings: &findings)
        if let hostContext {
            evaluateHost(
                hostContext,
                fingerprint: hostFingerprint,
                environmentCount: max(Int(experiment?.environmentCount ?? 1), 1),
                irreversibleProgram: irreversibleProgram,
                claims: &claims,
                findings: &findings
            )
        } else {
            claims.append(.init(
                identifier: "SC-HOST-001",
                title: "Host context is structurally valid",
                status: .notEvaluated,
                rationale: "No HostContextPack was supplied."
            ))
            limitations.append("Host, tissue, delivery, immune, metabolic, and containment context was not evaluated.")
        }

        let measuredOutputs: [String]
        if let experiment {
            measuredOutputs = evaluateExperiment(
                experiment,
                programPack: programPack,
                outputNames: outputNames,
                monitors: monitors,
                claims: &claims,
                findings: &findings
            )
        } else {
            measuredOutputs = []
            claims.append(.init(
                identifier: "SC-MEASURE-001",
                title: "Declared outputs are measured",
                status: .notEvaluated,
                rationale: "No ExperimentPack was supplied."
            ))
            limitations.append("Interventions, measurements, replicates, stop conditions, and checkpoints were not evaluated.")
        }

        let couplingCoverage: (required: [String], transactional: [String])
        if let coupling {
            couplingCoverage = evaluateCoupling(
                coupling,
                programPack: programPack,
                claims: &claims,
                findings: &findings
            )
        } else {
            couplingCoverage = ([], [])
            let f4 = programPack.header.fidelity == .tissue
            claims.append(.init(
                identifier: "SC-TXN-001",
                title: "Required external coupling is transactional",
                status: f4 ? .notEvaluated : .notApplicable,
                rationale: f4 ? "F4 was requested but no CouplingPack was supplied." : "No coupling contract was supplied."
            ))
            if f4 {
                add(&findings, "SC-F4-COUPLING", .error, "fidelity", "F4",
                    "Tissue fidelity is declared without an external coupling contract.",
                    "Supply a fingerprint-bound CouplingPack before interpreting F4 results.")
            }
        }

        let disposition: VivoSafetyDisposition = findings.contains { $0.severity == .blocking }
            ? .rejected
            : (findings.contains { $0.severity == .error }
                ? .conditionalComputationalUse
                : .computationalResearchOnly)
        return VivoSafetyCase(
            generatedAt: generatedAt,
            programFingerprint: programPack.header.contentFingerprint,
            programSourceFingerprint: programPack.header.sourceFingerprint,
            hostContextFingerprint: hostFingerprint,
            experimentFingerprint: experimentFingerprint,
            couplingFingerprint: couplingFingerprint,
            disposition: disposition,
            claims: claims.sorted { $0.identifier < $1.identifier },
            findings: findings.sorted {
                let lhs = severityRank($0.severity)
                let rhs = severityRank($1.severity)
                return lhs == rhs ? $0.identifier < $1.identifier : lhs > rhs
            },
            coverage: VivoSafetyCoverage(
                outputSpecies: outputNames,
                monitoredOutputSpecies: monitoredOutputs,
                measuredOutputSpecies: measuredOutputs,
                terminationMonitors: termination,
                requiredCouplingChannels: couplingCoverage.required,
                transactionalRequiredCouplingChannels: couplingCoverage.transactional,
                parameterEvidenceCounts: parameterCoverage
            ),
            limitations: limitations
        )
    }

    private func monitorReferences(
        monitors: [VivoProgramPack.MonitorMetadata],
        instructions: [VivoProgramPack.ExpressionInstructionMetadata],
        speciesCount: Int,
        findings: inout [VivoSafetyFinding]
    ) -> Set<UInt32> {
        var result = Set<UInt32>()
        for monitor in monitors {
            let start = Int(monitor.expressionOffset)
            let end = start.addingReportingOverflow(Int(monitor.expressionCount))
            guard !end.overflow, start >= 0, end.partialValue <= instructions.count else {
                add(&findings, "SC-MON-INVALID-\(monitor.index)", .blocking, "monitor-integrity", monitor.identifier,
                    "Monitor expression range is outside the ProgramPack expression table.",
                    "Reject and regenerate the ProgramPack with a compatible compiler.")
                continue
            }
            for instruction in instructions[start..<end.partialValue]
                where instruction.opcode == .loadSpecies && instruction.operand < UInt32(speciesCount) {
                result.insert(instruction.operand)
            }
        }
        return result
    }

    private func evaluateParameters(
        _ parameters: [VivoProgramPack.ParameterMetadata],
        claims: inout [VivoSafetyClaim],
        findings: inout [VivoSafetyFinding]
    ) -> [String: UInt32] {
        var counts: [String: UInt32] = [:]
        var unbounded: [String] = []
        var weakEvidence: [String] = []
        for parameter in parameters {
            counts[evidenceLabel(parameter.evidenceClass), default: 0] &+= 1
            let bounded = parameter.minimum.isFinite && parameter.maximum.isFinite &&
                parameter.minimum > -Double.greatestFiniteMagnitude &&
                parameter.maximum < Double.greatestFiniteMagnitude &&
                parameter.minimum <= parameter.value && parameter.value <= parameter.maximum
            if !bounded { unbounded.append(parameter.identifier) }
            if parameter.evidenceClass >= 4 { weakEvidence.append(parameter.identifier) }
        }
        claims.append(.init(
            identifier: "SC-PARAM-001",
            title: "Parameter uncertainty and evidence are explicit",
            status: unbounded.isEmpty && weakEvidence.isEmpty ? .satisfied : .conditional,
            rationale: "Parameter ranges and evidence classes were evaluated from the ProgramPack.",
            evidence: counts.keys.sorted().map { "\($0):\(counts[$0] ?? 0)" },
            assumptions: (unbounded + weakEvidence).sorted()
        ))
        if !unbounded.isEmpty {
            add(&findings, "SC-PARAM-BOUNDS", .error, "parameter-uncertainty", "parameters",
                "Parameters lack finite ordered ranges containing nominal values: \(unbounded.sorted().joined(separator: ", ")).",
                "Provide uncertainty intervals and propagate them through calibration or ensemble analysis.")
        }
        if !weakEvidence.isEmpty {
            add(&findings, "SC-PARAM-EVIDENCE", .warning, "parameter-evidence", "parameters",
                "Assumed or hypothetical parameters remain: \(weakEvidence.sorted().joined(separator: ", ")).",
                "Retain these classifications and replace them with stronger evidence before stronger claims.")
        }
        return counts
    }

    private func evaluateHost(
        _ host: VivoHostContextPack,
        fingerprint: VivoFingerprint?,
        environmentCount: Int,
        irreversibleProgram: Bool,
        claims: inout [VivoSafetyClaim],
        findings: inout [VivoSafetyFinding]
    ) {
        do {
            try host.validate(environmentCount: environmentCount)
            claims.append(.init(
                identifier: "SC-HOST-001",
                title: "Host context is structurally valid",
                status: .satisfied,
                rationale: "HostContextPack validation passed.",
                evidence: fingerprint.map { [$0.hex] } ?? []
            ))
        } catch {
            claims.append(.init(
                identifier: "SC-HOST-001",
                title: "Host context is structurally valid",
                status: .unsupported,
                rationale: String(describing: error)
            ))
            add(&findings, "SC-HOST-INVALID", .blocking, "host-context", "hostContext",
                String(describing: error), "Correct the HostContextPack before execution.")
        }

        let delivery = host.delivery
        let containment = delivery.externalShutdownSupported ||
            !delivery.molecularContainment.isEmpty ||
            (delivery.mode == .retrievableImplant && delivery.retrievalSupported)
        if delivery.isPersistent || irreversibleProgram {
            claims.append(.init(
                identifier: "SC-CONT-001",
                title: "Persistent or irreversible behavior has declared containment",
                status: containment ? .conditional : .unsupported,
                rationale: containment
                    ? "A containment path is declared; this software does not establish its reliability."
                    : "No independent containment path is declared.",
                evidence: delivery.physicalContainment + delivery.molecularContainment,
                assumptions: delivery.assumptions
            ))
            if !containment {
                add(&findings, "SC-CONT-MISSING", .blocking, "containment", delivery.mode.rawValue,
                    "Persistent or irreversible behavior lacks an independent containment path.",
                    "Provide retrieval, external shutdown, or molecular containment with quantified failure assumptions.")
            }
        } else {
            claims.append(.init(
                identifier: "SC-CONT-001",
                title: "Persistent or irreversible behavior has declared containment",
                status: .notApplicable,
                rationale: "Neither the program nor delivery context is classified as persistent or irreversible."
            ))
        }
        if delivery.isPersistent && delivery.locality == .systemic {
            add(&findings, "SC-CONT-SYSTEMIC", .error, "distribution", delivery.mode.rawValue,
                "Persistent delivery is declared systemic.",
                "Model distribution, persistence, shedding, off-target exposure, and independent shutdown.")
        }
        if host.uncertaintyTags.isEmpty {
            add(&findings, "SC-HOST-UNCERTAINTY", .warning, "uncertainty", "hostContext",
                "Host context declares no named uncertainty ranges.",
                "Attach uncertainty intervals to influential host, tissue, delivery, and transport assumptions.")
        }
    }

    private func evaluateExperiment(
        _ experiment: VivoExperimentPack,
        programPack: VivoProgramPack,
        outputNames: [String],
        monitors: [VivoProgramPack.MonitorMetadata],
        claims: inout [VivoSafetyClaim],
        findings: inout [VivoSafetyFinding]
    ) -> [String] {
        if experiment.programFingerprint != programPack.header.contentFingerprint {
            add(&findings, "SC-EXP-FINGERPRINT", .blocking, "artifact-identity", "experiment",
                "ExperimentPack references a different ProgramPack.",
                "Use an ExperimentPack bound to this ProgramPack fingerprint.")
        }
        let measured = Set(experiment.measurements.map { $0.species })
        let covered = outputNames.filter { measured.contains($0) }.sorted()
        let missing = outputNames.filter { !measured.contains($0) }.sorted()
        claims.append(.init(
            identifier: "SC-MEASURE-001",
            title: "Declared outputs are measured",
            status: outputNames.isEmpty ? .notApplicable : (missing.isEmpty ? .satisfied : .unsupported),
            rationale: "\(covered.count) of \(outputNames.count) outputs appear in the measurement plan.",
            evidence: covered,
            assumptions: missing
        ))
        for name in missing {
            add(&findings, "SC-MEASURE-MISSING-\(name)", .error, "measurement-coverage", name,
                "Declared output is absent from the experiment measurement plan.",
                "Add a bounded-cadence measurement or a justified validated downstream observable.")
        }
        let monitorNames = Set(monitors.map { $0.identifier })
        for condition in experiment.stopConditions
            where condition.required && !monitorNames.contains(condition.monitorIdentifier) {
            add(&findings, "SC-STOP-\(condition.identifier)", .blocking, "stop-condition", condition.identifier,
                "Required stop condition references unknown monitor \(condition.monitorIdentifier).",
                "Resolve the monitor against the same ProgramPack.")
        }
        if experiment.checkpointCadence == nil {
            add(&findings, "SC-CHECKPOINT", .warning, "recoverability", "experiment",
                "Experiment does not declare checkpoint cadence.",
                "Add bounded checkpoint cadence for long-running or expensive campaigns.")
        }
        return covered
    }

    private func evaluateCoupling(
        _ coupling: VivoCouplingPack,
        programPack: VivoProgramPack,
        claims: inout [VivoSafetyClaim],
        findings: inout [VivoSafetyFinding]
    ) -> (required: [String], transactional: [String]) {
        if coupling.programFingerprint != programPack.header.contentFingerprint {
            add(&findings, "SC-COUPLING-FINGERPRINT", .blocking, "artifact-identity", "coupling",
                "CouplingPack references a different ProgramPack.",
                "Use a CouplingPack compiled for this ProgramPack fingerprint.")
        }
        let required = coupling.channels.filter { $0.required }.map { $0.identifier }.sorted()
        let transactional = coupling.channels.filter { $0.required && $0.transactional }.map { $0.identifier }.sorted()
        let nonTransactional = coupling.channels.filter { $0.required && !$0.transactional }
        claims.append(.init(
            identifier: "SC-TXN-001",
            title: "Required external coupling is transactional",
            status: required.isEmpty ? .notApplicable : (nonTransactional.isEmpty ? .satisfied : .unsupported),
            rationale: "\(transactional.count) of \(required.count) required channels are transactional.",
            evidence: transactional,
            assumptions: nonTransactional.map { $0.identifier }
        ))
        for channel in nonTransactional {
            add(&findings, "SC-TXN-\(channel.identifier)", .blocking, "transactional-coupling", channel.identifier,
                "Required channel can cross the global transaction boundary independently.",
                "Mark it transactional and reject the global step when any required participant rejects.")
        }
        if coupling.policy.requireAllParticipants &&
            (coupling.policy.failureAction == .holdLastValue || coupling.policy.failureAction == .isolateSubsystem) {
            add(&findings, "SC-COUPLING-POLICY", .error, "coupling-failure-policy", "coupling.policy",
                "Failure policy can continue with stale or isolated state despite requiring all participants.",
                "Use rejectGlobalStep, substep, or an explicit shutdown response.")
        }
        for channel in coupling.channels where channel.required && channel.destinationBounds == nil {
            add(&findings, "SC-COUPLING-BOUND-\(channel.identifier)", .warning, "coupling-bound", channel.identifier,
                "Required channel has no explicit destination bounds.",
                "Declare destination bounds and enforce them before commit.")
        }
        return (required, transactional)
    }

    private func fingerprint<T: Codable & Sendable>(_ value: T?) throws -> VivoFingerprint? {
        guard let value else { return nil }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(value))
    }

    private func add(
        _ findings: inout [VivoSafetyFinding],
        _ id: String,
        _ severity: VivoSafetyFinding.Severity,
        _ category: String,
        _ subject: String,
        _ message: String,
        _ mitigation: String
    ) {
        findings.append(.init(
            identifier: id,
            severity: severity,
            category: category,
            subject: subject,
            message: message,
            requiredMitigation: mitigation
        ))
    }

    private func evidenceLabel(_ raw: UInt32) -> String {
        switch raw {
        case 0: "observed"
        case 1: "derived"
        case 2: "calibrated"
        case 3: "inferred"
        case 4: "assumed"
        case 5: "hypothetical"
        default: "unknown-\(raw)"
        }
    }

    private func severityRank(_ severity: VivoSafetyFinding.Severity) -> Int {
        switch severity {
        case .note: 0
        case .warning: 1
        case .error: 2
        case .blocking: 3
        }
    }
}
