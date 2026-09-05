import Foundation

public struct VivoSBMLQuantityMapping: Codable, Equatable, Sendable {
    public let unit: String
    public let scale: Double
    public let offset: Double
    public let minimum: Double
    public let maximum: Double

    public init(
        unit: String,
        scale: Double = 1,
        offset: Double = 0,
        minimum: Double,
        maximum: Double
    ) {
        self.unit = unit
        self.scale = scale
        self.offset = offset
        self.minimum = minimum
        self.maximum = maximum
    }

    public func transform(_ value: Double) throws -> Double {
        guard !unit.isEmpty,
              [scale, offset, minimum, maximum].allSatisfy(\.isFinite),
              scale != 0,
              minimum <= maximum else {
            throw VivoSBMLMappingError.invalidQuantityMapping(unit)
        }
        let transformed = value * scale + offset
        guard transformed.isFinite,
              transformed >= minimum,
              transformed <= maximum else {
            throw VivoSBMLMappingError.transformedValueOutsideBounds(
                unit: unit,
                value: transformed
            )
        }
        return transformed
    }
}

public struct VivoSBMLProgramMappingOptions: Codable, Equatable, Sendable {
    public var programName: String
    public var programVersion: String
    public var targetCellType: String
    public var targetTissue: String
    public var targetSpecies: String
    public var targetDiseaseState: String
    public var deliveryMode: String
    public var minimumFidelity: VivoFidelity
    public var maximumDurationSeconds: Double
    public var speciesMappings: [String: VivoSBMLQuantityMapping]
    public var parameterMappings: [String: VivoSBMLQuantityMapping]
    public var permitBoundarySpeciesAsInputs: Bool
    public var requireIntegerStoichiometry: Bool
    public var requireRecognizedKinetics: Bool

    public init(
        programName: String,
        programVersion: String = "0.1.0",
        targetCellType: String = "imported-sbml-cell",
        targetTissue: String = "imported-sbml-compartment-system",
        targetSpecies: String = "computational-host",
        targetDiseaseState: String = "not-declared",
        deliveryMode: String = "digital-model",
        minimumFidelity: VivoFidelity = .deterministic,
        maximumDurationSeconds: Double = 86_400,
        speciesMappings: [String: VivoSBMLQuantityMapping] = [:],
        parameterMappings: [String: VivoSBMLQuantityMapping] = [:],
        permitBoundarySpeciesAsInputs: Bool = true,
        requireIntegerStoichiometry: Bool = true,
        requireRecognizedKinetics: Bool = true
    ) {
        self.programName = programName
        self.programVersion = programVersion
        self.targetCellType = targetCellType
        self.targetTissue = targetTissue
        self.targetSpecies = targetSpecies
        self.targetDiseaseState = targetDiseaseState
        self.deliveryMode = deliveryMode
        self.minimumFidelity = minimumFidelity
        self.maximumDurationSeconds = maximumDurationSeconds
        self.speciesMappings = speciesMappings
        self.parameterMappings = parameterMappings
        self.permitBoundarySpeciesAsInputs = permitBoundarySpeciesAsInputs
        self.requireIntegerStoichiometry = requireIntegerStoichiometry
        self.requireRecognizedKinetics = requireRecognizedKinetics
    }

    public func validate() throws {
        guard Self.validIdentifier(programName),
              !programVersion.isEmpty,
              !targetCellType.isEmpty,
              !targetTissue.isEmpty,
              !targetSpecies.isEmpty,
              !deliveryMode.isEmpty,
              maximumDurationSeconds.isFinite,
              maximumDurationSeconds > 0 else {
            throw VivoSBMLMappingError.invalidOptions
        }
        for mapping in Array(speciesMappings.values) + Array(parameterMappings.values) {
            _ = try mapping.transform(mapping.minimum)
            _ = try mapping.transform(mapping.maximum)
        }
    }

    private static func validIdentifier(_ value: String) -> Bool {
        guard let first = value.unicodeScalars.first,
              CharacterSet.letters.union(CharacterSet(charactersIn: "_")).contains(first) else {
            return false
        }
        let allowed = CharacterSet.alphanumerics.union(CharacterSet(charactersIn: "_.-"))
        return value.unicodeScalars.dropFirst().allSatisfy(allowed.contains)
    }
}

public enum VivoSBMLMappingIssueSeverity: String, Codable, CaseIterable, Sendable {
    case note
    case warning
    case error
}

public struct VivoSBMLMappingIssue: Codable, Equatable, Sendable {
    public let severity: VivoSBMLMappingIssueSeverity
    public let code: String
    public let path: String
    public let message: String
}

public struct VivoSBMLProgramMappingResult: Sendable {
    public let programJSON: Data?
    public let issues: [VivoSBMLMappingIssue]
    public let sourceFingerprint: String
    public let normalizedModelFingerprint: String
    public let mappedSpeciesCount: Int
    public let mappedParameterCount: Int
    public let mappedReactionCount: Int
    public let unmappedReactionCount: Int

    public var succeeded: Bool {
        programJSON != nil && !issues.contains(where: { $0.severity == .error })
    }
}

public enum VivoSBMLMappingError: Error, LocalizedError, Sendable {
    case invalidOptions
    case invalidQuantityMapping(String)
    case transformedValueOutsideBounds(unit: String, value: Double)
    case invalidSourceModel
    case encoding(String)

    public var errorDescription: String? {
        switch self {
        case .invalidOptions: return "SBML mapping options are invalid"
        case .invalidQuantityMapping(let unit): return "SBML quantity mapping is invalid for unit \(unit)"
        case .transformedValueOutsideBounds(let unit, let value): return "mapped value \(value) is outside bounds for unit \(unit)"
        case .invalidSourceModel: return "SBML import result does not contain a valid normalized model"
        case .encoding(let reason): return "VivoProgram JSON encoding failed: \(reason)"
        }
    }
}

public struct VivoSBMLProgramMapper: Sendable {
    public init() {}

    public func map(
        _ imported: VivoSBMLImportResult,
        options: VivoSBMLProgramMappingOptions
    ) throws -> VivoSBMLProgramMappingResult {
        try options.validate()
        guard let model = imported.model else {
            throw VivoSBMLMappingError.invalidSourceModel
        }

        var issues: [VivoSBMLMappingIssue] = imported.issues.map {
            .init(
                severity: $0.severity == .error ? .error : ($0.severity == .warning ? .warning : .note),
                code: $0.code,
                path: $0.path,
                message: $0.message
            )
        }
        var inputObjects: [[String: Any]] = []
        var speciesObjects: [[String: Any]] = []
        var parameterObjects: [[String: Any]] = []
        var reactionObjects: [[String: Any]] = []
        var mappedSpecies: [String: MappedSpecies] = [:]
        var mappedParameters: [String: MappedParameter] = [:]

        for (index, species) in model.species.enumerated() {
            let path = "$.species[\(index)]"
            guard let mapping = resolveSpeciesMapping(species, options: options) else {
                issues.append(.init(
                    severity: .error,
                    code: "NVSMP001",
                    path: path,
                    message: "Species '\(species.id)' requires an explicit source-to-NumiVivo unit and bound mapping."
                ))
                continue
            }
            let sourceInitial = species.initialAmount ?? species.initialConcentration ?? 0
            let initial: Double
            do { initial = try mapping.transform(sourceInitial) }
            catch {
                issues.append(.init(
                    severity: .error,
                    code: "NVSMP002",
                    path: path,
                    message: error.localizedDescription
                ))
                continue
            }
            let isInput = species.boundaryCondition || species.constant
            if isInput && !options.permitBoundarySpeciesAsInputs {
                issues.append(.init(
                    severity: .error,
                    code: "NVSMP003",
                    path: path,
                    message: "Boundary or constant species '\(species.id)' cannot be mapped under the current policy."
                ))
                continue
            }
            mappedSpecies[species.id] = .init(
                id: species.id,
                unit: mapping.unit,
                isInput: isInput
            )
            let evidence: [String: Any] = [
                "class": "assumed",
                "source": model.sourceFingerprint,
                "context": "SBML import with explicit quantity mapping"
            ]
            if isInput {
                inputObjects.append([
                    "id": species.id,
                    "source": "experiment",
                    "unit": mapping.unit,
                    "default": initial,
                    "bounds": ["min": mapping.minimum, "max": mapping.maximum],
                    "evidence": evidence
                ])
            } else {
                speciesObjects.append([
                    "id": species.id,
                    "kind": mapping.unit == "count" ? "count" : "concentration",
                    "compartment": species.compartment,
                    "unit": mapping.unit,
                    "initial": initial,
                    "bounds": ["min": mapping.minimum, "max": mapping.maximum],
                    "conserved": false,
                    "externallyOwned": false,
                    "evidence": evidence
                ])
            }
        }

        for (index, parameter) in model.parameters.enumerated() {
            mapParameter(
                parameter,
                promotedID: parameter.id,
                path: "$.parameters[\(index)]",
                options: options,
                model: model,
                mappedParameters: &mappedParameters,
                parameterObjects: &parameterObjects,
                issues: &issues
            )
        }
        for (reactionIndex, reaction) in model.reactions.enumerated() {
            for (parameterIndex, parameter) in (reaction.kineticLaw?.localParameters ?? []).enumerated() {
                let promoted = Self.promotedParameterID(reaction: reaction.id, local: parameter.id)
                mapParameter(
                    parameter,
                    promotedID: promoted,
                    path: "$.reactions[\(reactionIndex)].kineticLaw.localParameters[\(parameterIndex)]",
                    options: options,
                    model: model,
                    mappedParameters: &mappedParameters,
                    parameterObjects: &parameterObjects,
                    issues: &issues,
                    mappingLookupID: "\(reaction.id).\(parameter.id)"
                )
            }
        }

        var unmappedReactions = 0
        for (reactionIndex, reaction) in model.reactions.enumerated() {
            let path = "$.reactions[\(reactionIndex)]"
            guard reaction.reactants.allSatisfy({ mappedSpecies[$0.species] != nil }),
                  reaction.products.allSatisfy({ mappedSpecies[$0.species] != nil }),
                  reaction.modifiers.allSatisfy({ mappedSpecies[$0.species] != nil }) else {
                issues.append(.init(
                    severity: .error,
                    code: "NVSMP010",
                    path: path,
                    message: "Reaction references at least one species that was not mapped."
                ))
                unmappedReactions += 1
                continue
            }
            guard let law = reaction.kineticLaw else {
                issues.append(.init(
                    severity: .error,
                    code: "NVSMP011",
                    path: path,
                    message: "Reaction has no kinetic law."
                ))
                unmappedReactions += 1
                continue
            }
            let expression: RateExpression?
            if let formula = law.formula, !formula.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
                var parser = FormulaParser(formula)
                expression = parser.parse()
            } else if let mathML = law.mathML {
                expression = MathMLExpressionParser.parse(mathML)
            } else {
                expression = nil
            }
            guard let expression else {
                issues.append(.init(
                    severity: .error,
                    code: "NVSMP012",
                    path: path,
                    message: "Kinetic law could not be parsed into the supported arithmetic subset."
                ))
                unmappedReactions += 1
                continue
            }

            let localMap = Dictionary(uniqueKeysWithValues: (law.localParameters).map {
                ($0.id, Self.promotedParameterID(reaction: reaction.id, local: $0.id))
            })
            let mappedLaw = matchRateLaw(
                expression,
                reaction: reaction,
                localParameters: localMap,
                globalParameters: mappedParameters,
                options: options
            )
            guard let mappedLaw else {
                issues.append(.init(
                    severity: options.requireRecognizedKinetics ? .error : .warning,
                    code: "NVSMP013",
                    path: path,
                    message: "Kinetic law is not an exact supported zero-order, mass-action, or reversible mass-action form."
                ))
                unmappedReactions += 1
                continue
            }

            var reactants: [[String: Any]] = []
            var products: [[String: Any]] = []
            var stoichiometryValid = true
            for reference in reaction.reactants {
                guard let coefficient = integerStoichiometry(reference.stoichiometry) else {
                    issues.append(.init(
                        severity: .error,
                        code: "NVSMP014",
                        path: path,
                        message: "Reactant stoichiometry for '\(reference.species)' is not a supported positive Int16."
                    ))
                    stoichiometryValid = false
                    continue
                }
                reactants.append(["species": reference.species, "coefficient": coefficient])
            }
            for reference in reaction.products {
                guard let coefficient = integerStoichiometry(reference.stoichiometry) else {
                    issues.append(.init(
                        severity: .error,
                        code: "NVSMP015",
                        path: path,
                        message: "Product stoichiometry for '\(reference.species)' is not a supported positive Int16."
                    ))
                    stoichiometryValid = false
                    continue
                }
                products.append(["species": reference.species, "coefficient": coefficient])
            }
            guard stoichiometryValid else {
                unmappedReactions += 1
                continue
            }
            reactionObjects.append([
                "id": reaction.id,
                "compartment": reaction.compartment ?? inferredCompartment(reaction, species: model.species),
                "reactants": reactants,
                "products": products,
                "rate": [
                    "law": mappedLaw.name,
                    "parameters": mappedLaw.parameters
                ],
                "critical": false
            ])
        }

        inputObjects.sort { ($0["id"] as? String ?? "") < ($1["id"] as? String ?? "") }
        speciesObjects.sort { ($0["id"] as? String ?? "") < ($1["id"] as? String ?? "") }
        parameterObjects.sort { ($0["id"] as? String ?? "") < ($1["id"] as? String ?? "") }
        reactionObjects.sort { ($0["id"] as? String ?? "") < ($1["id"] as? String ?? "") }

        let hasErrors = issues.contains(where: { $0.severity == .error })
        let programJSON: Data?
        if hasErrors {
            programJSON = nil
        } else {
            let program: [String: Any] = [
                "apiVersion": "numivivo.org/v1alpha1",
                "kind": "VivoProgram",
                "metadata": [
                    "name": options.programName,
                    "version": options.programVersion,
                    "description": "Conservative mass-action mapping from SBML source.",
                    "namespace": "urn:numivivo:sbml-import:\(model.fingerprint)",
                    "labels": [
                        "sbml-source": model.sourceFingerprint,
                        "sbml-normalized-model": model.fingerprint,
                        "mapping": "explicit-unit-and-supported-kinetics-only"
                    ]
                ],
                "spec": [
                    "target": [
                        "cellType": options.targetCellType,
                        "tissue": options.targetTissue,
                        "species": options.targetSpecies,
                        "developmentalStage": "not-declared",
                        "diseaseState": options.targetDiseaseState,
                        "deliveryMode": options.deliveryMode
                    ],
                    "minimumFidelity": Self.fidelityName(options.minimumFidelity),
                    "inputs": inputObjects,
                    "species": speciesObjects,
                    "state": [],
                    "parameters": parameterObjects,
                    "reactions": reactionObjects,
                    "rules": [],
                    "constraints": [],
                    "termination": [[
                        "id": "imported-model-maximum-duration",
                        "when": [
                            "gte": [
                                ["time": true],
                                ["literal": ["value": options.maximumDurationSeconds, "unit": "s"]]
                            ]
                        ],
                        "action": "permanent-shutdown",
                        "reason": "The bounded imported-model duration ended."
                    ]]
                ]
            ]
            do {
                programJSON = try JSONSerialization.data(
                    withJSONObject: program,
                    options: [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
                )
            } catch {
                throw VivoSBMLMappingError.encoding(error.localizedDescription)
            }
        }

        return .init(
            programJSON: programJSON,
            issues: issues,
            sourceFingerprint: model.sourceFingerprint,
            normalizedModelFingerprint: model.fingerprint,
            mappedSpeciesCount: mappedSpecies.count,
            mappedParameterCount: mappedParameters.count,
            mappedReactionCount: reactionObjects.count,
            unmappedReactionCount: unmappedReactions
        )
    }

    private struct MappedSpecies {
        let id: String
        let unit: String
        let isInput: Bool
    }

    private struct MappedParameter {
        let id: String
        let unit: String
        let value: Double
    }

    private struct MappedRateLaw {
        let name: String
        let parameters: [String]
    }

    private func resolveSpeciesMapping(
        _ species: VivoSBMLSpecies,
        options: VivoSBMLProgramMappingOptions
    ) -> VivoSBMLQuantityMapping? {
        if let explicit = options.speciesMappings[species.id] { return explicit }
        let unit = species.substanceUnits
        if species.hasOnlySubstanceUnits,
           unit == "item" || unit == "count" {
            return .init(unit: "count", minimum: 0, maximum: Double(UInt32.max))
        }
        if let unit, Self.recognizedNumiVivoUnits.contains(unit) {
            return .init(unit: unit, minimum: 0, maximum: Double.greatestFiniteMagnitude)
        }
        return nil
    }

    private func mapParameter(
        _ parameter: VivoSBMLParameter,
        promotedID: String,
        path: String,
        options: VivoSBMLProgramMappingOptions,
        model: VivoSBMLModel,
        mappedParameters: inout [String: MappedParameter],
        parameterObjects: inout [[String: Any]],
        issues: inout [VivoSBMLMappingIssue],
        mappingLookupID: String? = nil
    ) {
        let mapping = options.parameterMappings[mappingLookupID ?? parameter.id] ??
                      options.parameterMappings[parameter.id] ??
                      inferredParameterMapping(parameter)
        guard let mapping else {
            issues.append(.init(
                severity: .error,
                code: "NVSMP004",
                path: path,
                message: "Parameter '\(parameter.id)' requires an explicit source-to-NumiVivo unit and bound mapping."
            ))
            return
        }
        guard let sourceValue = parameter.value else {
            issues.append(.init(
                severity: .error,
                code: "NVSMP005",
                path: path,
                message: "Parameter '\(parameter.id)' has no finite value."
            ))
            return
        }
        let value: Double
        do { value = try mapping.transform(sourceValue) }
        catch {
            issues.append(.init(
                severity: .error,
                code: "NVSMP006",
                path: path,
                message: error.localizedDescription
            ))
            return
        }
        if mappedParameters[promotedID] != nil {
            issues.append(.init(
                severity: .error,
                code: "NVSMP007",
                path: path,
                message: "Mapped parameter identifier '\(promotedID)' is duplicated."
            ))
            return
        }
        mappedParameters[promotedID] = .init(id: promotedID, unit: mapping.unit, value: value)
        parameterObjects.append([
            "id": promotedID,
            "unit": mapping.unit,
            "value": value,
            "bounds": ["min": mapping.minimum, "max": mapping.maximum],
            "evidence": [
                "class": "assumed",
                "source": model.sourceFingerprint,
                "context": "SBML import with explicit quantity mapping"
            ]
        ])
    }

    private func inferredParameterMapping(
        _ parameter: VivoSBMLParameter
    ) -> VivoSBMLQuantityMapping? {
        guard let unit = parameter.units,
              Self.recognizedNumiVivoUnits.contains(unit),
              let value = parameter.value else { return nil }
        let magnitude = max(abs(value), 1e-12)
        return .init(
            unit: unit,
            minimum: value >= 0 ? 0 : -magnitude * 1000,
            maximum: magnitude * 1000
        )
    }

    private func matchRateLaw(
        _ expression: RateExpression,
        reaction: VivoSBMLReaction,
        localParameters: [String: String],
        globalParameters: [String: MappedParameter],
        options: VivoSBMLProgramMappingOptions
    ) -> MappedRateLaw? {
        let parameterIDs = Set(globalParameters.keys).union(localParameters.keys)
        func promoted(_ symbol: String) -> String? {
            if let local = localParameters[symbol] { return local }
            return globalParameters[symbol] != nil ? symbol : nil
        }
        let reactantPowers = stoichiometricPowers(reaction.reactants)
        let productPowers = stoichiometricPowers(reaction.products)

        if reaction.reversible,
           case .subtract(let forward, let reverse) = expression,
           let left = monomial(forward, parameterIDs: parameterIDs),
           let right = monomial(reverse, parameterIDs: parameterIDs),
           left.speciesPowers == reactantPowers,
           right.speciesPowers == productPowers,
           let forwardParameter = promoted(left.parameter),
           let reverseParameter = promoted(right.parameter) {
            return .init(
                name: "reversible-mass-action",
                parameters: [forwardParameter, reverseParameter]
            )
        }

        guard !reaction.reversible || reaction.products.isEmpty,
              let forward = monomial(expression, parameterIDs: parameterIDs),
              forward.speciesPowers == reactantPowers,
              let parameter = promoted(forward.parameter) else {
            return nil
        }
        return .init(
            name: reaction.reactants.isEmpty ? "zero-order" : "mass-action",
            parameters: [parameter]
        )
    }

    private func monomial(
        _ expression: RateExpression,
        parameterIDs: Set<String>
    ) -> (parameter: String, speciesPowers: [String: Int])? {
        var factors: [RateExpression] = []
        expression.flattenMultiplication(into: &factors)
        var parameter: String?
        var speciesPowers: [String: Int] = [:]
        for factor in factors {
            switch factor {
            case .number(let value):
                guard abs(value - 1) <= 1e-12 else { return nil }
            case .symbol(let symbol):
                if parameterIDs.contains(symbol) {
                    guard parameter == nil else { return nil }
                    parameter = symbol
                } else {
                    speciesPowers[symbol, default: 0] += 1
                }
            case .power(let base, let exponent):
                guard case .symbol(let symbol) = base,
                      case .number(let numericExponent) = exponent,
                      numericExponent.isFinite,
                      numericExponent >= 1,
                      numericExponent.rounded() == numericExponent,
                      numericExponent <= Double(Int.max),
                      !parameterIDs.contains(symbol) else {
                    return nil
                }
                speciesPowers[symbol, default: 0] += Int(numericExponent)
            default:
                return nil
            }
        }
        guard let parameter else { return nil }
        return (parameter, speciesPowers)
    }

    private func stoichiometricPowers(
        _ references: [VivoSBMLSpeciesReference]
    ) -> [String: Int] {
        var result: [String: Int] = [:]
        for reference in references {
            guard let value = integerStoichiometry(reference.stoichiometry) else { continue }
            result[reference.species, default: 0] += value
        }
        return result
    }

    private func integerStoichiometry(_ value: Double) -> Int? {
        guard value.isFinite,
              value > 0,
              value.rounded() == value,
              value <= Double(Int16.max) else { return nil }
        return Int(value)
    }

    private func inferredCompartment(
        _ reaction: VivoSBMLReaction,
        species: [VivoSBMLSpecies]
    ) -> String {
        let byID = Dictionary(uniqueKeysWithValues: species.map { ($0.id, $0.compartment) })
        let compartments = Set(
            (reaction.reactants.map(\.species) + reaction.products.map(\.species)).compactMap { byID[$0] }
        )
        return compartments.count == 1 ? compartments.first! : "cell"
    }

    private static func promotedParameterID(reaction: String, local: String) -> String {
        let base = reaction + "__" + local
        return base.map { character in
            character.isLetter || character.isNumber || character == "_" ||
            character == "-" || character == "." ? character : "_"
        }.reduce(into: "") { $0.append($1) }
    }

    private static func fidelityName(_ fidelity: VivoFidelity) -> String {
        switch fidelity {
        case .logic: return "F0"
        case .deterministic: return "F1"
        case .stochastic: return "F2"
        case .spatial: return "F3"
        case .tissue: return "F4"
        }
    }

    private static let recognizedNumiVivoUnits: Set<String> = [
        "1", "dimensionless", "normalized", "fraction", "%",
        "count", "molecule", "copy",
        "M", "mM", "uM", "µM", "nM", "pM", "fM",
        "s", "ms", "us", "µs", "min", "h", "d",
        "Hz", "1/s", "s^-1", "1/min", "min^-1", "1/h", "h^-1",
        "M/s", "mM/s", "uM/s", "µM/s", "nM/s", "nM/min",
        "count/s", "count/min", "ng/s", "ng/min", "ng/hour",
        "M^-1 s^-1", "mM^-1 s^-1", "uM^-1 s^-1", "nM^-1 s^-1"
    ]
}

private indirect enum RateExpression: Equatable {
    case symbol(String)
    case number(Double)
    case add(RateExpression, RateExpression)
    case subtract(RateExpression, RateExpression)
    case multiply(RateExpression, RateExpression)
    case divide(RateExpression, RateExpression)
    case power(RateExpression, RateExpression)

    func flattenMultiplication(into output: inout [RateExpression]) {
        if case .multiply(let left, let right) = self {
            left.flattenMultiplication(into: &output)
            right.flattenMultiplication(into: &output)
        } else {
            output.append(self)
        }
    }
}

private struct FormulaParser {
    private enum Token: Equatable {
        case identifier(String)
        case number(Double)
        case plus
        case minus
        case multiply
        case divide
        case power
        case leftParenthesis
        case rightParenthesis
        case comma
        case end
    }

    private let tokens: [Token]
    private var index = 0

    init(_ source: String) {
        self.tokens = Self.tokenize(source) ?? [.end]
    }

    mutating func parse() -> RateExpression? {
        guard let expression = parseExpression(), current == .end else { return nil }
        return expression
    }

    private mutating func parseExpression() -> RateExpression? {
        guard var value = parseTerm() else { return nil }
        while true {
            if consume(.plus) {
                guard let right = parseTerm() else { return nil }
                value = .add(value, right)
            } else if consume(.minus) {
                guard let right = parseTerm() else { return nil }
                value = .subtract(value, right)
            } else {
                return value
            }
        }
    }

    private mutating func parseTerm() -> RateExpression? {
        guard var value = parsePower() else { return nil }
        while true {
            if consume(.multiply) {
                guard let right = parsePower() else { return nil }
                value = .multiply(value, right)
            } else if consume(.divide) {
                guard let right = parsePower() else { return nil }
                value = .divide(value, right)
            } else {
                return value
            }
        }
    }

    private mutating func parsePower() -> RateExpression? {
        guard var value = parsePrimary() else { return nil }
        if consume(.power) {
            guard let exponent = parsePower() else { return nil }
            value = .power(value, exponent)
        }
        return value
    }

    private mutating func parsePrimary() -> RateExpression? {
        switch current {
        case .number(let value):
            advance()
            return .number(value)
        case .identifier(let identifier):
            advance()
            if identifier.lowercased() == "pow", consume(.leftParenthesis) {
                guard let base = parseExpression(),
                      consume(.comma),
                      let exponent = parseExpression(),
                      consume(.rightParenthesis) else { return nil }
                return .power(base, exponent)
            }
            return .symbol(identifier)
        case .leftParenthesis:
            advance()
            guard let value = parseExpression(), consume(.rightParenthesis) else { return nil }
            return value
        case .minus:
            advance()
            guard let value = parsePrimary() else { return nil }
            return .multiply(.number(-1), value)
        default:
            return nil
        }
    }

    private var current: Token { tokens[min(index, tokens.count - 1)] }

    @discardableResult
    private mutating func consume(_ token: Token) -> Bool {
        guard current == token else { return false }
        advance()
        return true
    }

    private mutating func advance() { index = min(index + 1, tokens.count - 1) }

    private static func tokenize(_ source: String) -> [Token]? {
        let scalars = Array(source.unicodeScalars)
        var result: [Token] = []
        var index = 0
        while index < scalars.count {
            let scalar = scalars[index]
            if CharacterSet.whitespacesAndNewlines.contains(scalar) {
                index += 1
                continue
            }
            switch scalar.value {
            case 43: result.append(.plus); index += 1
            case 45: result.append(.minus); index += 1
            case 42: result.append(.multiply); index += 1
            case 47: result.append(.divide); index += 1
            case 94: result.append(.power); index += 1
            case 40: result.append(.leftParenthesis); index += 1
            case 41: result.append(.rightParenthesis); index += 1
            case 44: result.append(.comma); index += 1
            default:
                if CharacterSet.decimalDigits.contains(scalar) || scalar.value == 46 {
                    let start = index
                    index += 1
                    while index < scalars.count,
                          CharacterSet(charactersIn: "0123456789.eE+-").contains(scalars[index]) {
                        if (scalars[index].value == 43 || scalars[index].value == 45),
                           scalars[index - 1].value != 101,
                           scalars[index - 1].value != 69 {
                            break
                        }
                        index += 1
                    }
                    guard let value = Double(String(String.UnicodeScalarView(scalars[start..<index]))), value.isFinite else {
                        return nil
                    }
                    result.append(.number(value))
                } else if CharacterSet.letters.union(CharacterSet(charactersIn: "_")).contains(scalar) {
                    let start = index
                    index += 1
                    let allowed = CharacterSet.alphanumerics.union(CharacterSet(charactersIn: "_.-"))
                    while index < scalars.count, allowed.contains(scalars[index]) { index += 1 }
                    result.append(.identifier(String(String.UnicodeScalarView(scalars[start..<index]))))
                } else {
                    return nil
                }
            }
        }
        result.append(.end)
        return result
    }
}

private final class MathMLExpressionParser: NSObject, XMLParserDelegate {
    private struct Frame {
        let name: String
        var text: String = ""
        var children: [RateExpression] = []
        var operatorName: String?
    }

    private var stack: [Frame] = []
    private var result: RateExpression?
    private var failed = false

    static func parse(_ source: String) -> RateExpression? {
        guard let data = source.data(using: .utf8) else { return nil }
        let delegate = MathMLExpressionParser()
        let parser = XMLParser(data: data)
        parser.delegate = delegate
        parser.shouldProcessNamespaces = true
        parser.shouldResolveExternalEntities = false
        return parser.parse() && !delegate.failed ? delegate.result : nil
    }

    func parser(
        _ parser: XMLParser,
        didStartElement elementName: String,
        namespaceURI: String?,
        qualifiedName qName: String?,
        attributes attributeDict: [String: String] = [:]
    ) {
        let name = localName(elementName)
        if let operation = Self.operation(name), !stack.isEmpty {
            stack[stack.count - 1].operatorName = operation
        } else {
            stack.append(.init(name: name))
        }
    }

    func parser(_ parser: XMLParser, foundCharacters string: String) {
        guard !stack.isEmpty else { return }
        stack[stack.count - 1].text += string
    }

    func parser(
        _ parser: XMLParser,
        didEndElement elementName: String,
        namespaceURI: String?,
        qualifiedName qName: String?
    ) {
        let name = localName(elementName)
        if Self.operation(name) != nil { return }
        guard !stack.isEmpty else { failed = true; return }
        let frame = stack.removeLast()
        guard frame.name == name else { failed = true; return }
        let expression: RateExpression?
        switch name {
        case "ci":
            let value = frame.text.trimmingCharacters(in: .whitespacesAndNewlines)
            expression = value.isEmpty ? nil : .symbol(value)
        case "cn":
            let value = Double(frame.text.trimmingCharacters(in: .whitespacesAndNewlines))
            expression = value.flatMap { $0.isFinite ? .number($0) : nil }
        case "apply":
            expression = build(operation: frame.operatorName, operands: frame.children)
        case "math":
            expression = frame.children.count == 1 ? frame.children[0] : nil
        default:
            expression = frame.children.count == 1 ? frame.children[0] : nil
        }
        guard let expression else { failed = true; return }
        if stack.isEmpty { result = expression }
        else { stack[stack.count - 1].children.append(expression) }
    }

    private func build(operation: String?, operands: [RateExpression]) -> RateExpression? {
        guard let operation else { return nil }
        switch operation {
        case "plus": return fold(operands, with: RateExpression.add)
        case "times": return fold(operands, with: RateExpression.multiply)
        case "minus":
            guard operands.count == 2 else { return nil }
            return .subtract(operands[0], operands[1])
        case "divide":
            guard operands.count == 2 else { return nil }
            return .divide(operands[0], operands[1])
        case "power":
            guard operands.count == 2 else { return nil }
            return .power(operands[0], operands[1])
        default: return nil
        }
    }

    private func fold(
        _ values: [RateExpression],
        with combine: (RateExpression, RateExpression) -> RateExpression
    ) -> RateExpression? {
        guard let first = values.first else { return nil }
        return values.dropFirst().reduce(first, combine)
    }

    private static func operation(_ name: String) -> String? {
        ["plus", "minus", "times", "divide", "power"].contains(name) ? name : nil
    }

    private func localName(_ value: String) -> String {
        value.split(separator: ":").last.map(String.init) ?? value
    }
}
