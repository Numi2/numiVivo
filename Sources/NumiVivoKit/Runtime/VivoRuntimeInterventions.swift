import Foundation

/// Accepted-boundary mutations execute synchronously on the transaction owner.
/// There is no suspension between reading and replacing state, no second clock,
/// and no facade-local lifecycle that can disagree with the GPU transaction.
public extension VivoTransactionalMolecularRuntime {
    func apply(_ context: PreparedVivoHostContext) throws {
        guard context.programFingerprint == pack.header.contentFingerprint,
              context.environmentCount == configuration.environmentCount,
              context.laneCount == configuration.laneCount,
              context.parameterValues.count == Int(pack.runtimeContract.parameterCount) * Int(configuration.environmentCount),
              context.transport.count == Int(pack.runtimeContract.speciesCount) else {
            throw VivoRuntimeError.invalidConfiguration("host context identity or layout differs from runtime")
        }
        // checkpoint rejects stopped, in-flight and prepared transactions.
        let before = try checkpoint()
        var state = try VivoLittleEndianFP32.decode(before.stateFP32LE)
        try applyBoundaryCoupling(context.initialCoupling, to: &state)
        let replacement = replacingBoundary(before, state: state, parameters: context.parameterValues, transport: context.transport)
        try restore(VivoMolecularResumeCheckpoint(configuration: configuration, state: replacement))
    }

    func apply(intervention: PreparedVivoIntervention.Operation) throws {
        switch intervention {
        case .coupling(let updates):
            let before = try checkpoint()
            var state = try VivoLittleEndianFP32.decode(before.stateFP32LE)
            try applyBoundaryCoupling(updates, to: &state)
            try restore(VivoMolecularResumeCheckpoint(configuration: configuration, state: replacingBoundary(before, state: state)))
        case .parameter(let index, let environments, let value):
            guard index < pack.runtimeContract.parameterCount, value.isFinite, !environments.isEmpty,
                  Set(environments).count == environments.count,
                  environments.allSatisfy({ $0 < configuration.environmentCount }) else {
                throw VivoRuntimeError.invalidConfiguration("parameter intervention index, environments or value")
            }
            let bounds = try pack.parameterMetadata()[Int(index)]
            guard Double(value) >= bounds.minimum, Double(value) <= bounds.maximum else {
                throw VivoRuntimeError.invalidConfiguration("parameter intervention violates declared bounds")
            }
            let before = try checkpoint()
            var parameters = try VivoLittleEndianFP32.decode(before.parametersFP32LE)
            for environment in environments {
                parameters[Int(index) * Int(configuration.environmentCount) + Int(environment)] = value
            }
            try restore(VivoMolecularResumeCheckpoint(configuration: configuration, state: replacingBoundary(before, parameters: parameters)))
        case .transport(let index, let value):
            guard index < pack.runtimeContract.speciesCount else { throw VivoRuntimeError.invalidConfiguration("transport intervention index") }
            let before = try checkpoint()
            var transport = try VivoTransportRecordLE.decode(before.transportRecordLE)
            transport[Int(index)] = value
            try setTransport(transport)
        case .reversibleShutdown(let reason):
            try stopReversibly(reason: reason)
        case .permanentShutdown(let reason):
            guard !reason.isEmpty else { throw VivoRuntimeError.invalidConfiguration("empty shutdown reason") }
            stopPermanently(reason: reason)
        }
    }

    private func applyBoundaryCoupling(_ updates: [VivoCouplingUpdate], to values: inout [Float]) throws {
        let metadata = try pack.speciesMetadata(), stride = Int(configuration.laneCount)
        guard values.count == metadata.count * stride else { throw VivoRuntimeError.invalidConfiguration("initial coupling state layout") }
        var destinations = Set<UInt64>()
        for update in updates {
            guard update.speciesIndex < UInt32(metadata.count), update.laneIndex < configuration.laneCount,
                  update.value.isFinite,
                  destinations.insert(UInt64(update.speciesIndex) << 32 | UInt64(update.laneIndex)).inserted else {
                throw VivoRuntimeError.invalidConfiguration("initial coupling index, value or duplicated destination")
            }
            let item = metadata[Int(update.speciesIndex)]
            guard item.isExternallyOwned || item.isInput else {
                throw VivoRuntimeError.invalidConfiguration("host coupling cannot write internally owned state")
            }
            let i = Int(update.speciesIndex) * stride + Int(update.laneIndex)
            switch update.mode {
            case .replace: values[i] = update.value
            case .add: values[i] += update.value
            case .relaxHalfway: values[i] += 0.5 * (update.value - values[i])
            }
            guard values[i].isFinite, values[i] >= item.minimum, values[i] <= item.maximum,
                  !item.isCountValued || (values[i] >= 0 && values[i] <= 16_777_216 && values[i].rounded() == values[i]) else {
                throw VivoRuntimeError.invalidConfiguration("host coupling violates state bounds or exact-count representation")
            }
        }
    }

    private func replacingBoundary(_ c: VivoMolecularCheckpoint, state: [Float]? = nil, parameters: [Float]? = nil,
                                   transport: [VivoSpeciesTransportABI]? = nil) -> VivoMolecularCheckpoint {
        .init(programFingerprint: c.programFingerprint, sourceProgramFingerprint: c.sourceProgramFingerprint,
              fidelity: c.fidelity, seed: c.seed, stepIndex: c.stepIndex, absoluteTimeSeconds: c.absoluteTimeSeconds,
              speciesCount: c.speciesCount, laneCount: c.laneCount, parameterCount: c.parameterCount,
              parameterEnvironmentCount: c.parameterEnvironmentCount, temporalStateCount: c.temporalStateCount,
              stateFP32LE: state.map(VivoLittleEndianFP32.encode) ?? c.stateFP32LE,
              parametersFP32LE: parameters.map(VivoLittleEndianFP32.encode) ?? c.parametersFP32LE,
              temporalStateFP32LE: c.temporalStateFP32LE,
              transportRecordLE: transport.map(VivoTransportRecordLE.encode) ?? c.transportRecordLE,
              velocityFP32LE: c.velocityFP32LE, volumeFractionFP32LE: c.volumeFractionFP32LE)
    }
}
