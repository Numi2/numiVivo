import CryptoKit
import Foundation

public struct VivoSBMLImportLimits: Sendable, Equatable {
    public var maximumBytes: Int
    public var maximumElements: Int
    public var maximumDepth: Int
    public var maximumAttributesPerElement: Int
    public var maximumTextBytes: Int
    public var maximumMathMLBytes: Int
    public var maximumCompartments: Int
    public var maximumSpecies: Int
    public var maximumParameters: Int
    public var maximumReactions: Int
    public var maximumReferencesPerReaction: Int

    public init(
        maximumBytes: Int = 64 * 1_024 * 1_024,
        maximumElements: Int = 2_000_000,
        maximumDepth: Int = 256,
        maximumAttributesPerElement: Int = 256,
        maximumTextBytes: Int = 16 * 1_024 * 1_024,
        maximumMathMLBytes: Int = 4 * 1_024 * 1_024,
        maximumCompartments: Int = 65_536,
        maximumSpecies: Int = 1_000_000,
        maximumParameters: Int = 1_000_000,
        maximumReactions: Int = 1_000_000,
        maximumReferencesPerReaction: Int = 65_536
    ) {
        self.maximumBytes = maximumBytes
        self.maximumElements = maximumElements
        self.maximumDepth = maximumDepth
        self.maximumAttributesPerElement = maximumAttributesPerElement
        self.maximumTextBytes = maximumTextBytes
        self.maximumMathMLBytes = maximumMathMLBytes
        self.maximumCompartments = maximumCompartments
        self.maximumSpecies = maximumSpecies
        self.maximumParameters = maximumParameters
        self.maximumReactions = maximumReactions
        self.maximumReferencesPerReaction = maximumReferencesPerReaction
    }

    public func validate() throws {
        guard maximumBytes > 0,
              maximumElements > 0,
              maximumDepth > 0,
              maximumAttributesPerElement > 0,
              maximumTextBytes > 0,
              maximumMathMLBytes > 0,
              maximumCompartments > 0,
              maximumSpecies > 0,
              maximumParameters > 0,
              maximumReactions > 0,
              maximumReferencesPerReaction > 0 else {
            throw VivoSBMLImportError.invalidLimits
        }
    }
}

public enum VivoSBMLIssueSeverity: String, Codable, CaseIterable, Sendable {
    case note
    case warning
    case error
}

public struct VivoSBMLImportIssue: Codable, Equatable, Sendable {
    public let severity: VivoSBMLIssueSeverity
    public let code: String
    public let path: String
    public let message: String
}

public struct VivoSBMLUnitComponent: Codable, Equatable, Sendable {
    public let kind: String
    public let exponent: Double
    public let scale: Int32
    public let multiplier: Double
    public let offset: Double
}

public struct VivoSBMLUnitDefinition: Codable, Equatable, Sendable {
    public let id: String
    public let name: String?
    public let units: [VivoSBMLUnitComponent]
}

public struct VivoSBMLCompartment: Codable, Equatable, Sendable {
    public let id: String
    public let name: String?
    public let spatialDimensions: Double
    public let size: Double?
    public let units: String?
    public let constant: Bool
    public let outside: String?
    public let sboTerm: String?
}

public struct VivoSBMLSpecies: Codable, Equatable, Sendable {
    public let id: String
    public let name: String?
    public let compartment: String
    public let initialAmount: Double?
    public let initialConcentration: Double?
    public let substanceUnits: String?
    public let hasOnlySubstanceUnits: Bool
    public let boundaryCondition: Bool
    public let constant: Bool
    public let conversionFactor: String?
    public let sboTerm: String?
}

public struct VivoSBMLParameter: Codable, Equatable, Sendable {
    public let id: String
    public let name: String?
    public let value: Double?
    public let units: String?
    public let constant: Bool
    public let sboTerm: String?
}

public struct VivoSBMLSpeciesReference: Codable, Equatable, Sendable {
    public let species: String
    public let stoichiometry: Double
    public let constant: Bool
    public let id: String?
    public let sboTerm: String?
}

public struct VivoSBMLModifierReference: Codable, Equatable, Sendable {
    public let species: String
    public let id: String?
    public let sboTerm: String?
}

public struct VivoSBMLKineticLaw: Codable, Equatable, Sendable {
    public let formula: String?
    public let mathML: String?
    public let substanceUnits: String?
    public let timeUnits: String?
    public let localParameters: [VivoSBMLParameter]
}

public struct VivoSBMLReaction: Codable, Equatable, Sendable {
    public let id: String
    public let name: String?
    public let reversible: Bool
    public let fast: Bool?
    public let compartment: String?
    public let reactants: [VivoSBMLSpeciesReference]
    public let products: [VivoSBMLSpeciesReference]
    public let modifiers: [VivoSBMLModifierReference]
    public let kineticLaw: VivoSBMLKineticLaw?
    public let sboTerm: String?
}

public struct VivoSBMLModel: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let level: UInt32
    public let version: UInt32
    public let id: String?
    public let name: String?
    public let substanceUnits: String?
    public let timeUnits: String?
    public let extentUnits: String?
    public let volumeUnits: String?
    public let areaUnits: String?
    public let lengthUnits: String?
    public let conversionFactor: String?
    public let unitDefinitions: [VivoSBMLUnitDefinition]
    public let compartments: [VivoSBMLCompartment]
    public let species: [VivoSBMLSpecies]
    public let parameters: [VivoSBMLParameter]
    public let reactions: [VivoSBMLReaction]
    public let sourceFingerprint: String
    public let fingerprint: String
}

public struct VivoSBMLImportResult: Codable, Equatable, Sendable {
    public let model: VivoSBMLModel?
    public let issues: [VivoSBMLImportIssue]
    public let sourceFingerprint: String
    public let sourceBytes: Int

    public var succeeded: Bool {
        model != nil && !issues.contains(where: { $0.severity == .error })
    }
}

public enum VivoSBMLImportError: Error, LocalizedError, Sendable {
    case invalidLimits
    case documentTooLarge(actual: Int, maximum: Int)
    case forbiddenMarkup(String)
    case parse(String)
    case structuralLimit(String)
    case invalidModel(String)
    case encoding(String)

    public var errorDescription: String? {
        switch self {
        case .invalidLimits: return "SBML import limits are invalid"
        case .documentTooLarge(let actual, let maximum): return "SBML document has \(actual) bytes; the limit is \(maximum)"
        case .forbiddenMarkup(let kind): return "SBML document contains forbidden \(kind) markup"
        case .parse(let reason): return "SBML parsing failed: \(reason)"
        case .structuralLimit(let reason): return "SBML document exceeds a structural limit: \(reason)"
        case .invalidModel(let reason): return "SBML model is invalid: \(reason)"
        case .encoding(let reason): return "SBML model fingerprint encoding failed: \(reason)"
        }
    }
}

public struct VivoSBMLImporter: Sendable {
    public let limits: VivoSBMLImportLimits

    public init(limits: VivoSBMLImportLimits = .init()) throws {
        try limits.validate()
        self.limits = limits
    }

    public func importModel(from data: Data) throws -> VivoSBMLImportResult {
        guard data.count <= limits.maximumBytes else {
            throw VivoSBMLImportError.documentTooLarge(
                actual: data.count,
                maximum: limits.maximumBytes
            )
        }
        if Self.containsASCII(data, token: "<!doctype") {
            throw VivoSBMLImportError.forbiddenMarkup("DOCTYPE")
        }
        if Self.containsASCII(data, token: "<!entity") {
            throw VivoSBMLImportError.forbiddenMarkup("ENTITY")
        }

        let sourceFingerprint = Self.hex(SHA256.hash(data: data))
        let delegate = VivoSBMLParserDelegate(
            limits: limits,
            sourceFingerprint: sourceFingerprint
        )
        let parser = XMLParser(data: data)
        parser.delegate = delegate
        parser.shouldProcessNamespaces = true
        parser.shouldReportNamespacePrefixes = false
        parser.shouldResolveExternalEntities = false
        let parsed = parser.parse()
        if let failure = delegate.failure {
            throw failure
        }
        if !parsed {
            let description = parser.parserError?.localizedDescription ?? "unknown XML parser error"
            throw VivoSBMLImportError.parse(description)
        }
        return try delegate.finish(sourceBytes: data.count)
    }

    private static func containsASCII(_ data: Data, token: String) -> Bool {
        let needle = Array(token.utf8)
        guard !needle.isEmpty, data.count >= needle.count else { return false }
        let bytes = Array(data)
        for start in 0...(bytes.count - needle.count) {
            var matches = true
            for offset in needle.indices {
                let byte = bytes[start + offset]
                let lower = byte >= 65 && byte <= 90 ? byte + 32 : byte
                if lower != needle[offset] {
                    matches = false
                    break
                }
            }
            if matches { return true }
        }
        return false
    }

    fileprivate static func hex<S: Sequence>(_ bytes: S) -> String where S.Element == UInt8 {
        bytes.map { String(format: "%02x", $0) }.joined()
    }
}

private final class VivoSBMLParserDelegate: NSObject, XMLParserDelegate {
    private enum ReferenceContext {
        case none
        case reactants
        case products
        case modifiers
    }

    private struct ModelBuilder {
        var level: UInt32?
        var version: UInt32?
        var id: String?
        var name: String?
        var substanceUnits: String?
        var timeUnits: String?
        var extentUnits: String?
        var volumeUnits: String?
        var areaUnits: String?
        var lengthUnits: String?
        var conversionFactor: String?
        var unitDefinitions: [VivoSBMLUnitDefinition] = []
        var compartments: [VivoSBMLCompartment] = []
        var species: [VivoSBMLSpecies] = []
        var parameters: [VivoSBMLParameter] = []
        var reactions: [VivoSBMLReaction] = []
    }

    private struct ReactionBuilder {
        var id = ""
        var name: String?
        var reversible = true
        var fast: Bool?
        var compartment: String?
        var reactants: [VivoSBMLSpeciesReference] = []
        var products: [VivoSBMLSpeciesReference] = []
        var modifiers: [VivoSBMLModifierReference] = []
        var kineticFormula: String?
        var kineticMathML: String?
        var kineticSubstanceUnits: String?
        var kineticTimeUnits: String?
        var localParameters: [VivoSBMLParameter] = []
        var sboTerm: String?

        func build() -> VivoSBMLReaction {
            let hasKinetic = kineticFormula != nil || kineticMathML != nil ||
                             !localParameters.isEmpty || kineticSubstanceUnits != nil ||
                             kineticTimeUnits != nil
            return .init(
                id: id,
                name: name,
                reversible: reversible,
                fast: fast,
                compartment: compartment,
                reactants: reactants,
                products: products,
                modifiers: modifiers,
                kineticLaw: hasKinetic ? .init(
                    formula: kineticFormula,
                    mathML: kineticMathML,
                    substanceUnits: kineticSubstanceUnits,
                    timeUnits: kineticTimeUnits,
                    localParameters: localParameters
                ) : nil,
                sboTerm: sboTerm
            )
        }
    }

    private struct UnitDefinitionBuilder {
        var id = ""
        var name: String?
        var units: [VivoSBMLUnitComponent] = []
    }

    let limits: VivoSBMLImportLimits
    let sourceFingerprint: String
    var failure: VivoSBMLImportError?

    private var model = ModelBuilder()
    private var reaction: ReactionBuilder?
    private var unitDefinition: UnitDefinitionBuilder?
    private var referenceContext = ReferenceContext.none
    private var depth = 0
    private var elementCount = 0
    private var textBytes = 0
    private var mathDepth = 0
    private var mathBuffer = ""
    private var skipDepth = 0
    private var issues: [VivoSBMLImportIssue] = []
    private var unsupportedFeatures: Set<String> = []
    private var sawSBML = false
    private var sawModel = false

    init(limits: VivoSBMLImportLimits, sourceFingerprint: String) {
        self.limits = limits
        self.sourceFingerprint = sourceFingerprint
    }

    func parser(
        _ parser: XMLParser,
        didStartElement elementName: String,
        namespaceURI: String?,
        qualifiedName qName: String?,
        attributes attributeDict: [String: String] = [:]
    ) {
        guard failure == nil else { parser.abortParsing(); return }
        depth += 1
        elementCount += 1
        guard depth <= limits.maximumDepth else {
            fail(.structuralLimit("element depth"), parser: parser)
            return
        }
        guard elementCount <= limits.maximumElements else {
            fail(.structuralLimit("element count"), parser: parser)
            return
        }
        guard attributeDict.count <= limits.maximumAttributesPerElement else {
            fail(.structuralLimit("attributes per element"), parser: parser)
            return
        }
        for (key, value) in attributeDict {
            guard key.utf8.count <= 4096,
                  value.utf8.count <= limits.maximumTextBytes else {
                fail(.structuralLimit("attribute bytes"), parser: parser)
                return
            }
        }

        let name = localName(elementName)
        if skipDepth > 0 {
            skipDepth += 1
            return
        }
        if name == "notes" || name == "annotation" {
            skipDepth = 1
            recordUnsupported(name, path: path(name), severity: .note)
            return
        }

        if mathDepth > 0 || name == "math" {
            mathDepth += 1
            appendMathStart(name: name, attributes: attributeDict, parser: parser)
            return
        }

        switch name {
        case "sbml":
            sawSBML = true
            model.level = uint32(attributeDict["level"])
            model.version = uint32(attributeDict["version"])
        case "model":
            sawModel = true
            model.id = attributeDict["id"]
            model.name = attributeDict["name"]
            model.substanceUnits = attributeDict["substanceUnits"]
            model.timeUnits = attributeDict["timeUnits"]
            model.extentUnits = attributeDict["extentUnits"]
            model.volumeUnits = attributeDict["volumeUnits"]
            model.areaUnits = attributeDict["areaUnits"]
            model.lengthUnits = attributeDict["lengthUnits"]
            model.conversionFactor = attributeDict["conversionFactor"]
        case "unitDefinition":
            unitDefinition = .init(
                id: attributeDict["id"] ?? "",
                name: attributeDict["name"],
                units: []
            )
        case "unit":
            guard unitDefinition != nil else { break }
            unitDefinition?.units.append(.init(
                kind: attributeDict["kind"] ?? "",
                exponent: finiteDouble(attributeDict["exponent"]) ?? 1,
                scale: int32(attributeDict["scale"]) ?? 0,
                multiplier: finiteDouble(attributeDict["multiplier"]) ?? 1,
                offset: finiteDouble(attributeDict["offset"]) ?? 0
            ))
        case "compartment":
            guard model.compartments.count < limits.maximumCompartments else {
                fail(.structuralLimit("compartment count"), parser: parser)
                return
            }
            model.compartments.append(.init(
                id: attributeDict["id"] ?? "",
                name: attributeDict["name"],
                spatialDimensions: finiteDouble(attributeDict["spatialDimensions"]) ?? 3,
                size: finiteDouble(attributeDict["size"]),
                units: attributeDict["units"],
                constant: boolean(attributeDict["constant"]) ?? true,
                outside: attributeDict["outside"],
                sboTerm: attributeDict["sboTerm"]
            ))
        case "species":
            guard model.species.count < limits.maximumSpecies else {
                fail(.structuralLimit("species count"), parser: parser)
                return
            }
            model.species.append(.init(
                id: attributeDict["id"] ?? "",
                name: attributeDict["name"],
                compartment: attributeDict["compartment"] ?? "",
                initialAmount: finiteDouble(attributeDict["initialAmount"]),
                initialConcentration: finiteDouble(attributeDict["initialConcentration"]),
                substanceUnits: attributeDict["substanceUnits"],
                hasOnlySubstanceUnits: boolean(attributeDict["hasOnlySubstanceUnits"]) ?? false,
                boundaryCondition: boolean(attributeDict["boundaryCondition"]) ?? false,
                constant: boolean(attributeDict["constant"]) ?? false,
                conversionFactor: attributeDict["conversionFactor"],
                sboTerm: attributeDict["sboTerm"]
            ))
        case "parameter":
            let parameter = VivoSBMLParameter(
                id: attributeDict["id"] ?? "",
                name: attributeDict["name"],
                value: finiteDouble(attributeDict["value"]),
                units: attributeDict["units"],
                constant: boolean(attributeDict["constant"]) ?? true,
                sboTerm: attributeDict["sboTerm"]
            )
            if reaction != nil {
                reaction?.localParameters.append(parameter)
            } else {
                guard model.parameters.count < limits.maximumParameters else {
                    fail(.structuralLimit("parameter count"), parser: parser)
                    return
                }
                model.parameters.append(parameter)
            }
        case "localParameter":
            guard reaction != nil else {
                issue(.error, "NVSBI001", path(name), "localParameter appears outside a kineticLaw")
                break
            }
            reaction?.localParameters.append(.init(
                id: attributeDict["id"] ?? "",
                name: attributeDict["name"],
                value: finiteDouble(attributeDict["value"]),
                units: attributeDict["units"],
                constant: true,
                sboTerm: attributeDict["sboTerm"]
            ))
        case "reaction":
            guard model.reactions.count < limits.maximumReactions else {
                fail(.structuralLimit("reaction count"), parser: parser)
                return
            }
            reaction = .init(
                id: attributeDict["id"] ?? "",
                name: attributeDict["name"],
                reversible: boolean(attributeDict["reversible"]) ?? true,
                fast: boolean(attributeDict["fast"]),
                compartment: attributeDict["compartment"],
                reactants: [],
                products: [],
                modifiers: [],
                kineticFormula: nil,
                kineticMathML: nil,
                kineticSubstanceUnits: nil,
                kineticTimeUnits: nil,
                localParameters: [],
                sboTerm: attributeDict["sboTerm"]
            )
        case "listOfReactants": referenceContext = .reactants
        case "listOfProducts": referenceContext = .products
        case "listOfModifiers": referenceContext = .modifiers
        case "speciesReference":
            guard var reaction else {
                issue(.error, "NVSBI002", path(name), "speciesReference appears outside a reaction")
                break
            }
            let reference = VivoSBMLSpeciesReference(
                species: attributeDict["species"] ?? "",
                stoichiometry: finiteDouble(attributeDict["stoichiometry"]) ?? 1,
                constant: boolean(attributeDict["constant"]) ?? true,
                id: attributeDict["id"],
                sboTerm: attributeDict["sboTerm"]
            )
            switch referenceContext {
            case .reactants:
                guard reaction.reactants.count < limits.maximumReferencesPerReaction else {
                    fail(.structuralLimit("reactants per reaction"), parser: parser)
                    return
                }
                reaction.reactants.append(reference)
            case .products:
                guard reaction.products.count < limits.maximumReferencesPerReaction else {
                    fail(.structuralLimit("products per reaction"), parser: parser)
                    return
                }
                reaction.products.append(reference)
            default:
                issue(.error, "NVSBI003", path(name), "speciesReference is not inside reactants or products")
            }
            self.reaction = reaction
        case "modifierSpeciesReference":
            guard var reaction else {
                issue(.error, "NVSBI004", path(name), "modifierSpeciesReference appears outside a reaction")
                break
            }
            guard reaction.modifiers.count < limits.maximumReferencesPerReaction else {
                fail(.structuralLimit("modifiers per reaction"), parser: parser)
                return
            }
            reaction.modifiers.append(.init(
                species: attributeDict["species"] ?? "",
                id: attributeDict["id"],
                sboTerm: attributeDict["sboTerm"]
            ))
            self.reaction = reaction
        case "kineticLaw":
            reaction?.kineticFormula = attributeDict["formula"]
            reaction?.kineticSubstanceUnits = attributeDict["substanceUnits"]
            reaction?.kineticTimeUnits = attributeDict["timeUnits"]
        case "assignmentRule", "rateRule", "algebraicRule", "initialAssignment",
             "event", "functionDefinition", "constraint", "compartmentType",
             "speciesType":
            recordUnsupported(name, path: path(name), severity: .warning)
        default:
            if name.hasPrefix("listOf") || name == "model" || name == "sbml" {
                break
            }
            if namespaceURI != nil,
               namespaceURI?.contains("sbml") == false,
               namespaceURI?.contains("MathML") == false {
                recordUnsupported("package:\(namespaceURI ?? "unknown")", path: path(name), severity: .warning)
            }
        }
    }

    func parser(
        _ parser: XMLParser,
        didEndElement elementName: String,
        namespaceURI: String?,
        qualifiedName qName: String?
    ) {
        guard failure == nil else { return }
        let name = localName(elementName)
        if skipDepth > 0 {
            skipDepth -= 1
            depth -= 1
            return
        }
        if mathDepth > 0 {
            appendMathEnd(name: name, parser: parser)
            mathDepth -= 1
            if mathDepth == 0 {
                reaction?.kineticMathML = mathBuffer
                mathBuffer = ""
            }
            depth -= 1
            return
        }

        switch name {
        case "unitDefinition":
            if let value = unitDefinition {
                model.unitDefinitions.append(.init(
                    id: value.id,
                    name: value.name,
                    units: value.units
                ))
            }
            unitDefinition = nil
        case "reaction":
            if let reaction { model.reactions.append(reaction.build()) }
            reaction = nil
            referenceContext = .none
        case "listOfReactants", "listOfProducts", "listOfModifiers":
            referenceContext = .none
        default: break
        }
        depth -= 1
    }

    func parser(_ parser: XMLParser, foundCharacters string: String) {
        guard failure == nil else { return }
        textBytes += string.utf8.count
        guard textBytes <= limits.maximumTextBytes else {
            fail(.structuralLimit("text bytes"), parser: parser)
            return
        }
        if mathDepth > 0 {
            appendMath(Self.escapeXML(string), parser: parser)
        }
    }

    func parser(_ parser: XMLParser, foundCDATA CDATABlock: Data) {
        guard failure == nil else { return }
        textBytes += CDATABlock.count
        guard textBytes <= limits.maximumTextBytes else {
            fail(.structuralLimit("CDATA bytes"), parser: parser)
            return
        }
        if mathDepth > 0 {
            appendMath(Self.escapeXML(String(decoding: CDATABlock, as: UTF8.self)), parser: parser)
        }
    }

    func parser(
        _ parser: XMLParser,
        resolveExternalEntityName name: String,
        systemID: String?
    ) -> Data? {
        issue(.error, "NVSBI005", "$", "External entity resolution was blocked: \(name)")
        return nil
    }

    func parser(_ parser: XMLParser, parseErrorOccurred parseError: Error) {
        if failure == nil {
            failure = .parse(parseError.localizedDescription)
        }
    }

    func finish(sourceBytes: Int) throws -> VivoSBMLImportResult {
        guard sawSBML else { throw VivoSBMLImportError.invalidModel("missing sbml root") }
        guard sawModel else { throw VivoSBMLImportError.invalidModel("missing model element") }
        guard let level = model.level, let version = model.version else {
            throw VivoSBMLImportError.invalidModel("missing level or version")
        }
        if level != 3 {
            issue(.warning, "NVSBI006", "$.sbml", "Importer preserves SBML Level \(level) but production mapping targets Level 3.")
        }
        validateReferences()

        let unsigned = VivoSBMLModel(
            schemaVersion: 1,
            level: level,
            version: version,
            id: model.id,
            name: model.name,
            substanceUnits: model.substanceUnits,
            timeUnits: model.timeUnits,
            extentUnits: model.extentUnits,
            volumeUnits: model.volumeUnits,
            areaUnits: model.areaUnits,
            lengthUnits: model.lengthUnits,
            conversionFactor: model.conversionFactor,
            unitDefinitions: model.unitDefinitions,
            compartments: model.compartments,
            species: model.species,
            parameters: model.parameters,
            reactions: model.reactions,
            sourceFingerprint: sourceFingerprint,
            fingerprint: ""
        )
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        let encoded: Data
        do { encoded = try encoder.encode(unsigned) }
        catch { throw VivoSBMLImportError.encoding(error.localizedDescription) }
        let fingerprint = VivoSBMLImporter.hex(SHA256.hash(data: encoded))
        let compiled = VivoSBMLModel(
            schemaVersion: unsigned.schemaVersion,
            level: unsigned.level,
            version: unsigned.version,
            id: unsigned.id,
            name: unsigned.name,
            substanceUnits: unsigned.substanceUnits,
            timeUnits: unsigned.timeUnits,
            extentUnits: unsigned.extentUnits,
            volumeUnits: unsigned.volumeUnits,
            areaUnits: unsigned.areaUnits,
            lengthUnits: unsigned.lengthUnits,
            conversionFactor: unsigned.conversionFactor,
            unitDefinitions: unsigned.unitDefinitions,
            compartments: unsigned.compartments,
            species: unsigned.species,
            parameters: unsigned.parameters,
            reactions: unsigned.reactions,
            sourceFingerprint: unsigned.sourceFingerprint,
            fingerprint: fingerprint
        )
        return .init(
            model: issues.contains(where: { $0.severity == .error }) ? nil : compiled,
            issues: issues,
            sourceFingerprint: sourceFingerprint,
            sourceBytes: sourceBytes
        )
    }

    private func validateReferences() {
        validateUnique(model.compartments.map(\.id), path: "$.compartments", code: "NVSBI010")
        validateUnique(model.species.map(\.id), path: "$.species", code: "NVSBI011")
        validateUnique(model.parameters.map(\.id), path: "$.parameters", code: "NVSBI012")
        validateUnique(model.reactions.map(\.id), path: "$.reactions", code: "NVSBI013")
        validateUnique(model.unitDefinitions.map(\.id), path: "$.unitDefinitions", code: "NVSBI014")

        let compartments = Set(model.compartments.map(\.id))
        let speciesIDs = Set(model.species.map(\.id))
        for (index, species) in model.species.enumerated() {
            if species.id.isEmpty {
                issue(.error, "NVSBI015", "$.species[\(index)].id", "Species identifier is empty.")
            }
            if !compartments.contains(species.compartment) {
                issue(.error, "NVSBI016", "$.species[\(index)].compartment", "Species references unknown compartment '\(species.compartment)'.")
            }
            if species.initialAmount != nil && species.initialConcentration != nil {
                issue(.error, "NVSBI017", "$.species[\(index)]", "Species declares both initialAmount and initialConcentration.")
            }
        }
        for (reactionIndex, reaction) in model.reactions.enumerated() {
            if reaction.id.isEmpty {
                issue(.error, "NVSBI018", "$.reactions[\(reactionIndex)].id", "Reaction identifier is empty.")
            }
            if let compartment = reaction.compartment, !compartments.contains(compartment) {
                issue(.error, "NVSBI019", "$.reactions[\(reactionIndex)].compartment", "Reaction references unknown compartment '\(compartment)'.")
            }
            for (referenceIndex, reference) in reaction.reactants.enumerated() {
                validateReference(reference.species, stoichiometry: reference.stoichiometry, known: speciesIDs, path: "$.reactions[\(reactionIndex)].reactants[\(referenceIndex)]")
            }
            for (referenceIndex, reference) in reaction.products.enumerated() {
                validateReference(reference.species, stoichiometry: reference.stoichiometry, known: speciesIDs, path: "$.reactions[\(reactionIndex)].products[\(referenceIndex)]")
            }
            for (referenceIndex, reference) in reaction.modifiers.enumerated() where !speciesIDs.contains(reference.species) {
                issue(.error, "NVSBI020", "$.reactions[\(reactionIndex)].modifiers[\(referenceIndex)]", "Modifier references unknown species '\(reference.species)'.")
            }
            if let law = reaction.kineticLaw {
                validateUnique(law.localParameters.map(\.id), path: "$.reactions[\(reactionIndex)].kineticLaw.localParameters", code: "NVSBI021")
                if law.formula == nil && law.mathML == nil {
                    issue(.warning, "NVSBI022", "$.reactions[\(reactionIndex)].kineticLaw", "Kinetic law contains no formula or MathML expression.")
                }
            } else {
                issue(.warning, "NVSBI023", "$.reactions[\(reactionIndex)]", "Reaction has no kinetic law.")
            }
        }
    }

    private func validateReference(
        _ species: String,
        stoichiometry: Double,
        known: Set<String>,
        path: String
    ) {
        if !known.contains(species) {
            issue(.error, "NVSBI024", path, "Reference names unknown species '\(species)'.")
        }
        if !stoichiometry.isFinite || stoichiometry <= 0 {
            issue(.error, "NVSBI025", path, "Stoichiometry must be finite and positive.")
        }
        if abs(stoichiometry.rounded() - stoichiometry) > 1e-9 {
            issue(.warning, "NVSBI026", path, "Non-integer stoichiometry requires an explicit NumiVivo mapping policy.")
        }
    }

    private func validateUnique(_ values: [String], path: String, code: String) {
        var seen: Set<String> = []
        for (index, value) in values.enumerated() {
            if value.isEmpty {
                issue(.error, code, "\(path)[\(index)]", "Required identifier is empty.")
            } else if !seen.insert(value).inserted {
                issue(.error, code, "\(path)[\(index)]", "Identifier '\(value)' is duplicated.")
            }
        }
    }

    private func appendMathStart(
        name: String,
        attributes: [String: String],
        parser: XMLParser
    ) {
        var fragment = "<\(name)"
        for key in attributes.keys.sorted() {
            fragment += " \(key)=\"\(Self.escapeXML(attributes[key] ?? ""))\""
        }
        fragment += ">"
        appendMath(fragment, parser: parser)
    }

    private func appendMathEnd(name: String, parser: XMLParser) {
        appendMath("</\(name)>", parser: parser)
    }

    private func appendMath(_ fragment: String, parser: XMLParser) {
        guard mathBuffer.utf8.count + fragment.utf8.count <= limits.maximumMathMLBytes else {
            fail(.structuralLimit("MathML bytes"), parser: parser)
            return
        }
        mathBuffer += fragment
    }

    private func recordUnsupported(
        _ feature: String,
        path: String,
        severity: VivoSBMLIssueSeverity
    ) {
        guard unsupportedFeatures.insert(feature).inserted else { return }
        issue(
            severity,
            "NVSBI100",
            path,
            "SBML feature '\(feature)' is preserved only as an import diagnostic and is not mapped into the current reaction-network model."
        )
    }

    private func issue(
        _ severity: VivoSBMLIssueSeverity,
        _ code: String,
        _ path: String,
        _ message: String
    ) {
        issues.append(.init(
            severity: severity,
            code: code,
            path: path,
            message: message
        ))
    }

    private func fail(_ error: VivoSBMLImportError, parser: XMLParser) {
        if failure == nil { failure = error }
        parser.abortParsing()
    }

    private func path(_ element: String) -> String {
        "$.xml.\(element)#\(elementCount)"
    }

    private func localName(_ value: String) -> String {
        value.split(separator: ":").last.map(String.init) ?? value
    }

    private func uint32(_ value: String?) -> UInt32? {
        guard let value, let parsed = UInt32(value) else { return nil }
        return parsed
    }

    private func int32(_ value: String?) -> Int32? {
        guard let value, let parsed = Int32(value) else { return nil }
        return parsed
    }

    private func finiteDouble(_ value: String?) -> Double? {
        guard let value, let parsed = Double(value), parsed.isFinite else { return nil }
        return parsed
    }

    private func boolean(_ value: String?) -> Bool? {
        guard let value else { return nil }
        switch value.lowercased() {
        case "true", "1": return true
        case "false", "0": return false
        default: return nil
        }
    }

    private static func escapeXML(_ value: String) -> String {
        value
            .replacingOccurrences(of: "&", with: "&amp;")
            .replacingOccurrences(of: "<", with: "&lt;")
            .replacingOccurrences(of: ">", with: "&gt;")
            .replacingOccurrences(of: "\"", with: "&quot;")
            .replacingOccurrences(of: "'", with: "&apos;")
    }
}
