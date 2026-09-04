import Foundation

public struct VivoSEDMLExportOptions: Sendable, Codable, Equatable {
    public var modelSource: String
    public var modelLanguage: String
    public var simulationIdentifier: String
    public var taskIdentifier: String
    public var reportIdentifier: String
    public var algorithmKISAOID: String

    public init(
        modelSource: String = "model.xml",
        modelLanguage: String = "urn:sedml:language:sbml.level-3.version-2",
        simulationIdentifier: String = "simulation",
        taskIdentifier: String = "task",
        reportIdentifier: String = "report",
        algorithmKISAOID: String = "KISAO:0000000"
    ) {
        self.modelSource = modelSource
        self.modelLanguage = modelLanguage
        self.simulationIdentifier = simulationIdentifier
        self.taskIdentifier = taskIdentifier
        self.reportIdentifier = reportIdentifier
        self.algorithmKISAOID = algorithmKISAOID
    }
}

public struct VivoSEDMLExport: Sendable {
    public let xml: Data
    public let fingerprint: VivoFingerprint
    public let diagnostics: [VivoStandardsDiagnostic]

    public init(
        xml: Data,
        fingerprint: VivoFingerprint,
        diagnostics: [VivoStandardsDiagnostic]
    ) {
        self.xml = xml
        self.fingerprint = fingerprint
        self.diagnostics = diagnostics
    }
}

public struct VivoSEDMLExporter: Sendable {
    public init() {}

    public func export(
        experiment: PreparedVivoExperiment,
        programPack: VivoProgramPack,
        identifiers: VivoIdentifierMap,
        options: VivoSEDMLExportOptions = .init()
    ) throws -> VivoSEDMLExport {
        guard experiment.programFingerprint == programPack.header.contentFingerprint else {
            throw VivoArtifactValidationError.incompatible("SED-ML experiment and ProgramPack fingerprints differ")
        }
        guard experiment.durationSeconds.isFinite,
              experiment.durationSeconds > 0,
              experiment.preferredTimeStepSeconds.isFinite,
              experiment.preferredTimeStepSeconds > 0 else {
            throw VivoArtifactValidationError.invalid("SED-ML export requires finite positive duration and preferred time step")
        }

        let species = try programPack.speciesMetadata()
        let nominalPoints = ceil(experiment.durationSeconds / experiment.preferredTimeStepSeconds)
        guard nominalPoints.isFinite, nominalPoints > 0, nominalPoints <= Double(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("SED-ML uniform time-course point count is outside UInt32")
        }
        let points = UInt32(nominalPoints)
        let simulationID = VivoXMLCodec.identifier(options.simulationIdentifier, fallback: "simulation")
        let taskID = VivoXMLCodec.identifier(options.taskIdentifier, fallback: "task")
        let reportID = VivoXMLCodec.identifier(options.reportIdentifier, fallback: "report")
        var diagnostics: [VivoStandardsDiagnostic] = []

        struct GeneratedMeasurement {
            let measurement: PreparedVivoMeasurement
            let speciesID: String
            let dataGeneratorID: String
            let variableID: String
        }
        var generated: [GeneratedMeasurement] = []
        generated.reserveCapacity(experiment.measurements.count)
        for (index, measurement) in experiment.measurements.enumerated() {
            guard Int(measurement.speciesIndex) < species.count else {
                throw VivoArtifactValidationError.invalid("SED-ML measurement \(measurement.identifier) references an invalid species index")
            }
            let speciesMetadata = species[Int(measurement.speciesIndex)]
            guard let speciesID = identifiers.species[speciesMetadata.identifier] else {
                throw VivoArtifactValidationError.unresolved("SED-ML has no SBML identifier for \(speciesMetadata.identifier)")
            }
            generated.append(.init(
                measurement: measurement,
                speciesID: speciesID,
                dataGeneratorID: VivoXMLCodec.identifier("dg_\(measurement.identifier)_\(index)", fallback: "dg_\(index)"),
                variableID: VivoXMLCodec.identifier("variable_\(measurement.identifier)_\(index)", fallback: "variable_\(index)")
            ))
            if measurement.lanes.count != 1 || String(describing: measurement.aggregation) != "direct" {
                diagnostics.append(.init(
                    severity: .warning,
                    code: "NVSEDML001",
                    subject: measurement.identifier,
                    message: "Lane selection and NumiVivo aggregation are retained in annotations; the SED-ML variable targets the underlying SBML species."
                ))
            }
        }

        var lines: [String] = []
        lines.reserveCapacity(80 + generated.count * 14)
        lines.append(#"<?xml version="1.0" encoding="UTF-8"?>"#)
        lines.append(#"<sedML xmlns="http://sed-ml.org/sed-ml/level1/version4" xmlns:numivivo="https://numivivo.org/ns#" xmlns:sbml="http://www.sbml.org/sbml/level3/version2/core" level="1" version="4">"#)
        lines.append("  <annotation>")
        lines.append(#"    <numivivo:experiment programFingerprint="\#(programPack.header.contentFingerprint.hex)" experimentFingerprint="\#(experiment.fingerprint.hex)" fidelity="\#(experiment.fidelity.label)" environmentCount="\#(experiment.environmentCount)" replicateCount="\#(experiment.replicateCount)" minimumTimeStep="\#(VivoXMLCodec.finite(experiment.minimumTimeStepSeconds))" maximumTimeStep="\#(VivoXMLCodec.finite(experiment.maximumTimeStepSeconds))"/>"#)
        lines.append("  </annotation>")
        lines.append("  <listOfModels>")
        lines.append(#"    <model id="model" name="NumiVivo \#(VivoXMLCodec.escapedAttribute(identifiers.model))" language="\#(VivoXMLCodec.escapedAttribute(options.modelLanguage))" source="\#(VivoXMLCodec.escapedAttribute(options.modelSource))"/>"#)
        lines.append("  </listOfModels>")
        lines.append("  <listOfSimulations>")
        lines.append(#"    <uniformTimeCourse id="\#(simulationID)" initialTime="0" outputStartTime="0" outputEndTime="\#(VivoXMLCodec.finite(experiment.durationSeconds))" numberOfPoints="\#(points)">"#)
        lines.append(#"      <algorithm kisaoID="\#(VivoXMLCodec.escapedAttribute(options.algorithmKISAOID))">"#)
        lines.append("        <annotation>")
        lines.append(#"          <numivivo:runtime algorithm="transactional-multifidelity" fidelity="\#(experiment.fidelity.label)" preferredTimeStep="\#(VivoXMLCodec.finite(experiment.preferredTimeStepSeconds))" adaptiveSubsteps="true"/>"#)
        lines.append("        </annotation>")
        lines.append("      </algorithm>")
        lines.append("    </uniformTimeCourse>")
        lines.append("  </listOfSimulations>")
        lines.append("  <listOfTasks>")
        lines.append(#"    <task id="\#(taskID)" modelReference="model" simulationReference="\#(simulationID)"/>"#)
        lines.append("  </listOfTasks>")

        if !generated.isEmpty {
            lines.append("  <listOfDataGenerators>")
            for item in generated {
                let measurement = item.measurement
                lines.append(#"    <dataGenerator id="\#(item.dataGeneratorID)" name="\#(VivoXMLCodec.escapedAttribute(measurement.identifier))">"#)
                lines.append("      <listOfVariables>")
                lines.append(#"        <variable id="\#(item.variableID)" taskReference="\#(taskID)" target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='\#(VivoXMLCodec.escapedAttribute(item.speciesID))']"/>"#)
                lines.append("      </listOfVariables>")
                lines.append(#"      <math xmlns="http://www.w3.org/1998/Math/MathML"><ci>\#(item.variableID)</ci></math>"#)
                lines.append("      <annotation>")
                lines.append(#"        <numivivo:measurement outputUnit="\#(VivoXMLCodec.escapedAttribute(measurement.outputUnit))" cadenceSeconds="\#(VivoXMLCodec.finite(measurement.cadenceSeconds))" startSeconds="\#(VivoXMLCodec.finite(measurement.startSeconds))" endSeconds="\#(VivoXMLCodec.finite(measurement.endSeconds))" lanes="\#(measurement.lanes.map(String.init).joined(separator: ","))" aggregation="\#(VivoXMLCodec.escapedAttribute(String(describing: measurement.aggregation)))" storage="\#(VivoXMLCodec.escapedAttribute(measurement.storage.rawValue))"/>"#)
                lines.append("      </annotation>")
                lines.append("    </dataGenerator>")
            }
            lines.append("  </listOfDataGenerators>")
            lines.append("  <listOfOutputs>")
            lines.append(#"    <report id="\#(reportID)" name="NumiVivo measurements">"#)
            lines.append("      <listOfDataSets>")
            for (index, item) in generated.enumerated() {
                let dataSetID = VivoXMLCodec.identifier("dataset_\(item.measurement.identifier)_\(index)", fallback: "dataset_\(index)")
                lines.append(#"        <dataSet id="\#(dataSetID)" label="\#(VivoXMLCodec.escapedAttribute(item.measurement.identifier))" dataReference="\#(item.dataGeneratorID)"/>"#)
            }
            lines.append("      </listOfDataSets>")
            lines.append("    </report>")
            lines.append("  </listOfOutputs>")
        }

        if !experiment.interventions.isEmpty {
            diagnostics.append(.init(
                severity: .warning,
                code: "NVSEDML002",
                subject: experiment.fingerprint.hex,
                message: "Timed NumiVivo interventions are not lowered into SED-ML model changes in this export; the authoritative ExperimentPack remains required for execution."
            ))
        }
        if !experiment.stopConditions.isEmpty {
            diagnostics.append(.init(
                severity: .warning,
                code: "NVSEDML003",
                subject: experiment.fingerprint.hex,
                message: "Transactional monitor stop responses are not expressible by a uniformTimeCourse task and remain authoritative in the ExperimentPack."
            ))
        }
        if experiment.replicateCount > 1 {
            diagnostics.append(.init(
                severity: .note,
                code: "NVSEDML004",
                subject: experiment.fingerprint.hex,
                message: "Replicate seeds and environment batching are recorded in annotation rather than expanded into repeatedTask elements."
            ))
        }
        diagnostics.append(.init(
            severity: .note,
            code: "NVSEDML005",
            subject: experiment.fingerprint.hex,
            message: "The generic KiSAO identifier is intentional: NumiVivo combines deterministic, stochastic, spatial, temporal, and transactional algorithms under one fidelity contract."
        ))

        lines.append("</sedML>")
        let text = lines.joined(separator: "\n") + "\n"
        guard let data = text.data(using: .utf8) else {
            throw VivoArtifactValidationError.invalid("SED-ML serialization did not produce UTF-8")
        }
        return VivoSEDMLExport(
            xml: data,
            fingerprint: try VivoCanonicalJSON.fingerprint(data),
            diagnostics: diagnostics
        )
    }
}
