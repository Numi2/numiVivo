import Foundation

public extension VivoProgramPack {
    enum RateLaw: UInt32, Codable, Sendable, CaseIterable {
        case zeroOrder = 0
        case massAction = 1
        case hillActivation = 2
        case hillRepression = 3
        case michaelisMenten = 4
        case reversibleMassAction = 5
        case passiveTransport = 6
        case saturableTransport = 7
        case degradation = 8
        case customBytecode = 255
    }

    enum ExpressionOpcode: UInt16, Codable, Sendable {
        case pushConstant = 0
        case loadSpecies = 1
        case loadParameter = 2
        case logicalNot = 3
        case logicalAnd = 4
        case logicalOr = 5
        case greater = 6
        case greaterEqual = 7
        case less = 8
        case lessEqual = 9
        case equal = 10
        case notEqual = 11
        case add = 12
        case subtract = 13
        case multiply = 14
        case divide = 15
        case minimum = 16
        case maximum = 17
        case clamp = 18
        case sustained = 19
        case within = 20
        case risingEdge = 21
        case fallingEdge = 22
        case end = 255
    }

    enum ActionKind: UInt32, Codable, Sendable {
        case setOutput = 0
        case addOutput = 1
        case express = 2
        case suppress = 3
        case degrade = 4
        case setState = 5
        case incrementState = 6
        case emitEvent = 7
        case requestDifferentiation = 8
        case requestMigration = 9
        case reversibleShutdown = 10
        case permanentShutdown = 11
    }

    enum MonitorResponse: UInt32, Codable, Sendable {
        case record = 0
        case clamp = 1
        case rejectStep = 2
        case substep = 3
        case reversibleShutdown = 4
        case permanentShutdown = 5
    }

    struct StoichiometryMetadata: Codable, Sendable, Equatable {
        public let speciesIndex: UInt32
        public let speciesIdentifier: String
        public let coefficient: Int16
        public let isProduct: Bool
    }

    struct ReactionMetadata: Codable, Sendable, Equatable {
        public let index: UInt32
        public let identifier: String
        public let compartment: String
        public let reactants: [StoichiometryMetadata]
        public let products: [StoichiometryMetadata]
        public let parameterIndices: [UInt32]
        public let parameterIdentifiers: [String]
        public let rateLaw: RateLaw
        public let flags: UInt32
        public let expressionOffset: UInt32
        public let expressionCount: UInt32
        public let gateExpressionOffset: UInt32?
        public let delaySeconds: Float
        public let characteristicRate: Float
        public let cohortIndex: UInt32

        public var isCritical: Bool { flags & (1 << 0) != 0 }
        public var hasGate: Bool { flags & (1 << 1) != 0 }
        public var isDelayed: Bool { flags & (1 << 2) != 0 }
        public var isStochasticEligible: Bool { flags & (1 << 3) != 0 }
        public var isSpatial: Bool { flags & (1 << 4) != 0 }
    }

    struct ExpressionInstructionMetadata: Codable, Sendable, Equatable {
        public let index: UInt32
        public let opcode: ExpressionOpcode
        public let flags: UInt16
        public let operand: UInt32
        public let immediate: Float
        public let auxiliary: UInt32
    }

    struct ActionMetadata: Codable, Sendable, Equatable {
        public let index: UInt32
        public let targetIndex: UInt32
        public let expressionOffset: UInt32
        public let expressionCount: UInt32
        public let kind: ActionKind
        public let constantValue: Float
        public let maximumRate: Float
        public let unit: String
        public let flags: UInt32

        public var targetIsStringOffset: Bool { flags & 1 != 0 }
    }

    struct RuleMetadata: Codable, Sendable, Equatable {
        public let index: UInt32
        public let identifier: String
        public let conditionOffset: UInt32
        public let conditionCount: UInt32
        public let actionOffset: UInt32
        public let actionCount: UInt32
        public let priority: Int32
        public let refractorySeconds: Float
        public let temporalStateOffset: UInt32?
    }

    struct MonitorMetadata: Codable, Sendable, Equatable {
        public let index: UInt32
        public let identifier: String
        public let expressionOffset: UInt32
        public let expressionCount: UInt32
        public let message: String
        public let severity: UInt32
        public let response: MonitorResponse
        public let temporalStateOffset: UInt32?
        public let flags: UInt32

        public var isTerminationRule: Bool { flags & 1 != 0 }
    }

    struct CohortMetadata: Codable, Sendable, Equatable {
        public let index: UInt32
        public let reactionOffset: UInt32
        public let reactionCount: UInt32
        public let rateLaw: RateLaw
        public let flags: UInt32
        public let maximumStableStep: Float
        public let stiffnessEstimate: Float
        public let preferredThreads: UInt32
    }

    struct SpeciesIncidenceMetadata: Codable, Sendable, Equatable {
        public let reactionIndex: UInt32
        public let netCoefficient: Int16
    }

    func stoichiometryMetadata() throws -> [StoichiometryMetadata] {
        let section = try section(.stoichiometry)
        guard section.stride == 8 else {
            throw VivoProgramPackError.invalidSection("stoichiometry stride must be 8 bytes")
        }
        let species = try speciesMetadata()
        let reader = VivoPackRecordReader(data: data)
        var result: [StoichiometryMetadata] = []
        result.reserveCapacity(Int(section.count))
        for index in 0..<Int(section.count) {
            let base = Int(section.offset) + index * Int(section.stride)
            let speciesIndex = try reader.u32(base)
            guard speciesIndex < species.count else {
                throw VivoProgramPackError.invalidSection("stoichiometry species index \(speciesIndex) is out of bounds")
            }
            let role = try reader.u16(base + 6)
            guard role <= 1 else {
                throw VivoProgramPackError.invalidSection("stoichiometry role must be reactant or product")
            }
            result.append(.init(
                speciesIndex: UInt32(speciesIndex),
                speciesIdentifier: species[speciesIndex].identifier,
                coefficient: try reader.i16(base + 4),
                isProduct: role == 1
            ))
        }
        return result
    }

    func reactionParameterIndexTable() throws -> [UInt32] {
        let section = try section(.reactionParameterIndices)
        guard section.stride == 4 else {
            throw VivoProgramPackError.invalidSection("reaction parameter index stride must be 4 bytes")
        }
        let reader = VivoPackRecordReader(data: data)
        return try (0..<Int(section.count)).map {
            try reader.u32(Int(section.offset) + $0 * 4)
        }
    }

    func reactionMetadata() throws -> [ReactionMetadata] {
        let section = try section(.reactions)
        guard section.stride == 64 else {
            throw VivoProgramPackError.invalidSection("reaction stride must be 64 bytes")
        }
        let species = try speciesMetadata()
        let parameters = try parameterMetadata()
        let stoichiometry = try stoichiometryMetadata()
        let parameterIndices = try reactionParameterIndexTable()
        let strings = try section(.strings)
        let reader = VivoPackRecordReader(data: data)
        var result: [ReactionMetadata] = []
        result.reserveCapacity(Int(section.count))

        for index in 0..<Int(section.count) {
            let base = Int(section.offset) + index * Int(section.stride)
            let reactantOffset = try reader.u32(base + 8)
            let productOffset = try reader.u32(base + 12)
            let reactantCount = try reader.u32(base + 16)
            let productCount = try reader.u32(base + 20)
            let parameterOffset = try reader.u32(base + 24)
            let parameterCount = try reader.u32(base + 28)
            let expressionOffset = try reader.u32(base + 32)
            let expressionCount = try reader.u32(base + 36)
            let rawRateLaw = try reader.u32(base + 40)
            guard let rateLaw = RateLaw(rawValue: rawRateLaw) else {
                throw VivoProgramPackError.invalidSection("reaction \(index) has unknown rate law \(rawRateLaw)")
            }
            let flags = try reader.u32(base + 44)
            let cohortIndex = try reader.u32(base + 56)
            let reserved = try reader.u32(base + 60)

            let reactantRange = try checkedRange(offset: reactantOffset, count: reactantCount, limit: stoichiometry.count, label: "reaction reactants")
            let productRange = try checkedRange(offset: productOffset, count: productCount, limit: stoichiometry.count, label: "reaction products")
            let parameterRange = try checkedRange(offset: parameterOffset, count: parameterCount, limit: parameterIndices.count, label: "reaction parameters")
            let reactionParameterIndices = Array(parameterIndices[parameterRange])
            for parameterIndex in reactionParameterIndices where parameterIndex >= parameters.count {
                throw VivoProgramPackError.invalidSection("reaction \(index) parameter index \(parameterIndex) is out of bounds")
            }
            let gateOffset: UInt32? = flags & (1 << 1) != 0 ? reserved : nil
            if expressionCount > 0 {
                _ = try checkedRange(offset: expressionOffset, count: expressionCount, limit: Int(try section(.expressions).count), label: "reaction expression")
            }
            if let gateOffset, gateOffset >= (try section(.expressions).count) {
                throw VivoProgramPackError.invalidSection("reaction \(index) gate offset is out of bounds")
            }

            result.append(.init(
                index: UInt32(index),
                identifier: try readString(offset: try reader.u32(base), strings: strings, reader: reader),
                compartment: try readString(offset: try reader.u32(base + 4), strings: strings, reader: reader),
                reactants: Array(stoichiometry[reactantRange]),
                products: Array(stoichiometry[productRange]),
                parameterIndices: reactionParameterIndices,
                parameterIdentifiers: reactionParameterIndices.map { parameters[$0].identifier },
                rateLaw: rateLaw,
                flags: flags,
                expressionOffset: expressionOffset,
                expressionCount: expressionCount,
                gateExpressionOffset: gateOffset,
                delaySeconds: try reader.f32(base + 48),
                characteristicRate: try reader.f32(base + 52),
                cohortIndex: cohortIndex
            ))
        }
        _ = species
        return result
    }

    func expressionInstructions() throws -> [ExpressionInstructionMetadata] {
        let section = try section(.expressions)
        guard section.stride == 16 else {
            throw VivoProgramPackError.invalidSection("expression instruction stride must be 16 bytes")
        }
        let reader = VivoPackRecordReader(data: data)
        var result: [ExpressionInstructionMetadata] = []
        result.reserveCapacity(Int(section.count))
        for index in 0..<Int(section.count) {
            let base = Int(section.offset) + index * Int(section.stride)
            let rawOpcode = try reader.u16(base)
            guard let opcode = ExpressionOpcode(rawValue: rawOpcode) else {
                throw VivoProgramPackError.invalidSection("expression instruction \(index) has unknown opcode \(rawOpcode)")
            }
            result.append(.init(
                index: UInt32(index),
                opcode: opcode,
                flags: try reader.u16(base + 2),
                operand: try reader.u32(base + 4),
                immediate: try reader.f32(base + 8),
                auxiliary: try reader.u32(base + 12)
            ))
        }
        return result
    }

    func actionMetadata() throws -> [ActionMetadata] {
        let section = try section(.actions)
        guard section.stride == 32 else {
            throw VivoProgramPackError.invalidSection("action stride must be 32 bytes")
        }
        let strings = try section(.strings)
        let reader = VivoPackRecordReader(data: data)
        var result: [ActionMetadata] = []
        result.reserveCapacity(Int(section.count))
        for index in 0..<Int(section.count) {
            let base = Int(section.offset) + index * Int(section.stride)
            let rawKind = try reader.u32(base + 12)
            guard let kind = ActionKind(rawValue: rawKind) else {
                throw VivoProgramPackError.invalidSection("action \(index) has unknown kind \(rawKind)")
            }
            result.append(.init(
                index: UInt32(index),
                targetIndex: try reader.u32(base),
                expressionOffset: try reader.u32(base + 4),
                expressionCount: try reader.u32(base + 8),
                kind: kind,
                constantValue: try reader.f32(base + 16),
                maximumRate: try reader.f32(base + 20),
                unit: try readString(offset: try reader.u32(base + 24), strings: strings, reader: reader),
                flags: try reader.u32(base + 28)
            ))
        }
        return result
    }

    func ruleMetadata() throws -> [RuleMetadata] {
        let section = try section(.rules)
        guard section.stride == 32 else {
            throw VivoProgramPackError.invalidSection("rule stride must be 32 bytes")
        }
        let strings = try section(.strings)
        let reader = VivoPackRecordReader(data: data)
        var result: [RuleMetadata] = []
        result.reserveCapacity(Int(section.count))
        for index in 0..<Int(section.count) {
            let base = Int(section.offset) + index * Int(section.stride)
            let temporal = try reader.u32(base + 28)
            result.append(.init(
                index: UInt32(index),
                identifier: try readString(offset: try reader.u32(base), strings: strings, reader: reader),
                conditionOffset: try reader.u32(base + 4),
                conditionCount: try reader.u32(base + 8),
                actionOffset: try reader.u32(base + 12),
                actionCount: try reader.u32(base + 16),
                priority: try reader.i32(base + 20),
                refractorySeconds: try reader.f32(base + 24),
                temporalStateOffset: temporal == UInt32.max ? nil : temporal
            ))
        }
        return result
    }

    func monitorMetadata() throws -> [MonitorMetadata] {
        let section = try section(.monitors)
        guard section.stride == 32 else {
            throw VivoProgramPackError.invalidSection("monitor stride must be 32 bytes")
        }
        let strings = try section(.strings)
        let reader = VivoPackRecordReader(data: data)
        var result: [MonitorMetadata] = []
        result.reserveCapacity(Int(section.count))
        for index in 0..<Int(section.count) {
            let base = Int(section.offset) + index * Int(section.stride)
            let rawResponse = try reader.u32(base + 20)
            guard let response = MonitorResponse(rawValue: rawResponse) else {
                throw VivoProgramPackError.invalidSection("monitor \(index) has unknown response \(rawResponse)")
            }
            let temporal = try reader.u32(base + 24)
            result.append(.init(
                index: UInt32(index),
                identifier: try readString(offset: try reader.u32(base), strings: strings, reader: reader),
                expressionOffset: try reader.u32(base + 4),
                expressionCount: try reader.u32(base + 8),
                message: try readString(offset: try reader.u32(base + 12), strings: strings, reader: reader),
                severity: try reader.u32(base + 16),
                response: response,
                temporalStateOffset: temporal == UInt32.max ? nil : temporal,
                flags: try reader.u32(base + 28)
            ))
        }
        return result
    }

    func cohortMetadata() throws -> [CohortMetadata] {
        let section = try section(.cohorts)
        guard section.stride == 32 else {
            throw VivoProgramPackError.invalidSection("cohort stride must be 32 bytes")
        }
        let reader = VivoPackRecordReader(data: data)
        var result: [CohortMetadata] = []
        result.reserveCapacity(Int(section.count))
        for index in 0..<Int(section.count) {
            let base = Int(section.offset) + index * Int(section.stride)
            let rawRateLaw = try reader.u32(base + 8)
            guard let rateLaw = RateLaw(rawValue: rawRateLaw) else {
                throw VivoProgramPackError.invalidSection("cohort \(index) has unknown rate law \(rawRateLaw)")
            }
            result.append(.init(
                index: UInt32(index),
                reactionOffset: try reader.u32(base),
                reactionCount: try reader.u32(base + 4),
                rateLaw: rateLaw,
                flags: try reader.u32(base + 12),
                maximumStableStep: try reader.f32(base + 16),
                stiffnessEstimate: try reader.f32(base + 20),
                preferredThreads: try reader.u32(base + 24)
            ))
        }
        return result
    }

    func speciesIncidenceMetadata() throws -> [[SpeciesIncidenceMetadata]] {
        let offsetsSection = try section(.speciesIncidenceOffsets)
        let incidenceSection = try section(.speciesIncidence)
        guard offsetsSection.stride == 4, incidenceSection.stride == 16 else {
            throw VivoProgramPackError.invalidSection("species incidence section strides are invalid")
        }
        guard offsetsSection.count == runtimeContract.speciesCount + 1 else {
            throw VivoProgramPackError.invalidSection("species incidence offsets must contain speciesCount + 1 entries")
        }
        let reader = VivoPackRecordReader(data: data)
        let offsets = try (0..<Int(offsetsSection.count)).map {
            try reader.u32(Int(offsetsSection.offset) + $0 * 4)
        }
        guard offsets.first == 0, offsets.last == incidenceSection.count else {
            throw VivoProgramPackError.invalidSection("species incidence offsets do not span the incidence table")
        }
        var result: [[SpeciesIncidenceMetadata]] = []
        result.reserveCapacity(Int(runtimeContract.speciesCount))
        for speciesIndex in 0..<Int(runtimeContract.speciesCount) {
            let begin = offsets[speciesIndex]
            let end = offsets[speciesIndex + 1]
            guard begin <= end, end <= incidenceSection.count else {
                throw VivoProgramPackError.invalidSection("species incidence offsets are not monotonic")
            }
            var links: [SpeciesIncidenceMetadata] = []
            links.reserveCapacity(Int(end - begin))
            for offset in begin..<end {
                let base = Int(incidenceSection.offset) + Int(offset) * Int(incidenceSection.stride)
                let reaction = try reader.u32(base)
                guard reaction < runtimeContract.reactionCount else {
                    throw VivoProgramPackError.invalidSection("species incidence reaction index is out of bounds")
                }
                links.append(.init(reactionIndex: reaction, netCoefficient: try reader.i16(base + 4)))
            }
            result.append(links)
        }
        return result
    }

    private func readString(
        offset: UInt32,
        strings: Section,
        reader: VivoPackRecordReader
    ) throws -> String {
        guard UInt64(offset) < strings.size else {
            throw VivoProgramPackError.invalidSection("string offset \(offset) is outside the string table")
        }
        var cursor = Int(strings.offset) + Int(offset)
        let limit = Int(strings.offset + strings.size)
        var bytes: [UInt8] = []
        while cursor < limit {
            let byte = try reader.u8(cursor)
            if byte == 0 { return String(decoding: bytes, as: UTF8.self) }
            bytes.append(byte)
            cursor += 1
        }
        throw VivoProgramPackError.invalidSection("unterminated string at offset \(offset)")
    }

    private func checkedRange(
        offset: UInt32,
        count: UInt32,
        limit: Int,
        label: String
    ) throws -> Range<Int> {
        let start = Int(offset)
        let result = start.addingReportingOverflow(Int(count))
        guard !result.overflow, start >= 0, result.partialValue <= limit else {
            throw VivoProgramPackError.invalidSection("\(label) range is out of bounds")
        }
        return start..<result.partialValue
    }
}

private struct VivoPackRecordReader: Sendable {
    let data: Data

    func u8(_ offset: Int) throws -> UInt8 {
        guard offset >= 0, offset < data.count else {
            throw VivoProgramPackError.truncated("record byte at offset \(offset)")
        }
        return data[data.index(data.startIndex, offsetBy: offset)]
    }

    func u16(_ offset: Int) throws -> UInt16 {
        UInt16(try u8(offset)) | UInt16(try u8(offset + 1)) << 8
    }

    func i16(_ offset: Int) throws -> Int16 {
        Int16(bitPattern: try u16(offset))
    }

    func u32(_ offset: Int) throws -> UInt32 {
        UInt32(try u8(offset)) |
        UInt32(try u8(offset + 1)) << 8 |
        UInt32(try u8(offset + 2)) << 16 |
        UInt32(try u8(offset + 3)) << 24
    }

    func i32(_ offset: Int) throws -> Int32 {
        Int32(bitPattern: try u32(offset))
    }

    func f32(_ offset: Int) throws -> Float {
        Float(bitPattern: try u32(offset))
    }
}
