import Foundation

/// Concrete ProgramPack execution contract, checked before GPU allocation.
/// Unrepresented state classes are rejected rather than silently discarded.
public enum VivoProgramExecutionContract {
    public static func validate(pack: VivoProgramPack, configuration: VivoRuntimeConfiguration) throws {
        try configuration.validate(for: pack)
        let inspection = VivoNativeCompilerBridge.inspectProgramPack(pack.data, verifySectionHashes: true)
        guard inspection.invocation.succeeded else {
            throw VivoRuntimeError.packError("ProgramPack failed native integrity/semantic inspection")
        }
        let strides: [(VivoProgramPack.SectionKind, UInt32)] = [
            (.species, 32), (.parameters, 48), (.reactionParameterIndices, 4),
            (.stoichiometry, 8), (.reactions, 64), (.expressions, 16),
            (.actions, 32), (.rules, 32), (.monitors, 32),
            (.speciesIncidenceOffsets, 4), (.speciesIncidence, 16)
        ]
        for (kind, expected) in strides {
            guard try pack.section(kind).stride == expected else {
                throw VivoRuntimeError.packError("unsupported \(kind) stride; expected \(expected)")
            }
        }
        let contract = pack.runtimeContract
        let offsetCount = try pack.section(.speciesIncidenceOffsets).count
        guard contract.speciesCount > 0, contract.maximumExpressionStack <= 256,
              UInt64(offsetCount) == UInt64(contract.speciesCount) + 1 else {
            throw VivoRuntimeError.packError("invalid species/incidence/stack shape")
        }
        for count in [contract.speciesCount, contract.reactionCount, contract.ruleCount,
                      contract.monitorCount, contract.temporalStateCount] {
            guard UInt64(count) * UInt64(configuration.laneCount) <= UInt64(UInt32.max) else {
                throw VivoRuntimeError.invalidConfiguration("runtime dispatch exceeds UInt32 addressing")
            }
        }
        let actionCount = try pack.section(.actions).count
        let possibleEventsPerLane = UInt64(actionCount) + UInt64(contract.monitorCount)
        let possibleEvents = possibleEventsPerLane.multipliedReportingOverflow(by: UInt64(configuration.laneCount))
        guard !possibleEvents.overflow, possibleEvents.partialValue <= UInt64(UInt32.max) else {
            throw VivoRuntimeError.invalidConfiguration("worst-case event counter would overflow UInt32")
        }
        if let grid = configuration.spatialGrid {
            guard grid.width <= UInt32(Int32.max), grid.height <= UInt32(Int32.max), grid.depth <= UInt32(Int32.max) else {
                throw VivoRuntimeError.invalidConfiguration("spatial coordinates exceed the signed neighbor-index ABI")
            }
        }
        let species = try pack.speciesMetadata()
        if configuration.fidelity.rawValue >= VivoFidelity.spatial.rawValue,
           species.contains(where: { $0.isCountValued }) {
            throw VivoRuntimeError.invalidConfiguration(
                "spatial count transport requires stochastic hopping; use concentration state for this finite-volume backend"
            )
        }
        let reactions = try pack.sectionData(.reactions)
        let expressions = try pack.sectionData(.expressions)
        let rules = try pack.sectionData(.rules)
        let stoichiometry = try pack.sectionData(.stoichiometry)
        func word(_ bytes: Data, _ offset: Int) -> UInt32 {
            bytes.withUnsafeBytes { UInt32(littleEndian: $0.loadUnaligned(fromByteOffset: offset, as: UInt32.self)) }
        }
        func checkPureExpression(start: UInt32, count: UInt32) throws {
            let total = expressions.count / 16
            guard Int(start) < total else { throw VivoRuntimeError.packError("expression offset out of range") }
            let remaining = total - Int(start)
            let length = count == 0 ? min(remaining, 4096) : Int(count)
            guard length <= remaining else { throw VivoRuntimeError.packError("expression range out of bounds") }
            var ended = false
            for i in 0..<length {
                let opcode = word(expressions, (Int(start) + i) * 16) & 0xffff
                if opcode == 255 { ended = true; break }
                guard !(19...22).contains(opcode) else {
                    throw VivoRuntimeError.invalidConfiguration(
                        "temporal operators in reaction rates/gates require a once-per-step evaluation stage; only rule/monitor temporal operators are supported"
                    )
                }
            }
            guard ended else { throw VivoRuntimeError.packError("unterminated expression") }
        }
        for index in 0..<Int(contract.reactionCount) {
            let base = index * 64
            let flags = word(reactions, base + 44)
            let delay = Float(bitPattern: word(reactions, base + 48))
            guard flags & 4 == 0, delay == 0 else {
                throw VivoRuntimeError.invalidConfiguration("delayed reactions require a queue-owning backend")
            }
            let expressionCount = word(reactions, base + 36)
            if expressionCount > 0 { try checkPureExpression(start: word(reactions, base + 32), count: expressionCount) }
            if flags & 2 != 0 { try checkPureExpression(start: word(reactions, base + 60), count: 0) }
            if configuration.fidelity == .stochastic {
                guard word(reactions, base + 40) != 255 else {
                    throw VivoRuntimeError.invalidConfiguration("custom bytecode has no declared stochastic propensity contract")
                }
                for (offsetField, countField) in [(8, 16), (12, 20)] {
                    let offset = word(reactions, base + offsetField)
                    let count = word(reactions, base + countField)
                    for term in 0..<count {
                        let speciesIndex = word(stoichiometry, Int(UInt64(offset) + UInt64(term)) * 8)
                        guard Int(speciesIndex) < species.count,
                              species[Int(speciesIndex)].isExternallyOwned || species[Int(speciesIndex)].isCountValued else {
                            throw VivoRuntimeError.invalidConfiguration(
                                "ProgramPack F2 requires count-valued reaction state; use an explicit hybrid plan for continuous components"
                            )
                        }
                    }
                }
            }
        }
        for index in 0..<Int(contract.ruleCount) {
            let refractory = Float(bitPattern: word(rules, index * 32 + 24))
            guard refractory == 0 else {
                throw VivoRuntimeError.invalidConfiguration("nonzero rule refractory time requires separately checkpointed refractory state")
            }
        }
    }
}
