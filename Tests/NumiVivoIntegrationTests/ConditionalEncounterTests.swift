import Foundation
import Testing
@testable import NumiVivoKit

/// Host tests of the actual encounter admission/probability owner. Synthetic
/// component descriptors do not claim a qualified TST source. Full positive
/// chemistry is exercised by the separately retained native whole-route test.
@Suite(.serialized) struct ConditionalEncounterTests {
    private func components(_ count: Int = 2) throws -> [VivoConditionalEncounterComponent] {
        try (0..<count).map { index in
            try .init(index: index, atomIdentifiers: ["atom-\(index)"],
                  qualifiedPointFingerprint: VivoCanonicalJSON.fingerprint(Data("host component \(index)".utf8)))
        }
    }

    private func request(_ components: [VivoConditionalEncounterComponent], tagged: Int = 1,
                         values: [Double] = [3], unit: VivoConditionalEncounterReservoirUnit = .pascal,
                         times: [Double] = [0, 0.125, 1], sensitivities: Bool = true,
                         temperature: Double = 298.15, endpoint: String = "reactants",
                         numerics: VivoKineticSensitivityConfiguration = .init()) -> VivoConditionalEncounterRequest {
        let others = components.filter { $0.index != tagged }
        return .init(identifier: "explicit maintained bath host conformance", reactantEndpointIdentifier: endpoint,
            temperatureK: temperature, taggedComponent: components[tagged],
            reservoirs: zip(others, values).map { .init(component: $0, value: $1, unit: unit) },
            observationTimesSeconds: times, includeReservoirSensitivities: sensitivities, numerics: numerics)
    }

    private func plan(_ request: VivoConditionalEncounterRequest, components: [VivoConditionalEncounterComponent],
                      logRate: Double = log(2), standardState: VivoThermoStandardState = .idealGas(pressurePa: 101325),
                      unitsOverride: String? = nil, molecularityOverride: Int? = nil) throws -> VivoConditionalEncounterPlan {
        let value: Double, standardUnits: String, rateUnits: String
        switch standardState {
        case .idealGas(let pressure):
            value = pressure; standardUnits = "Pa"; rateUnits = "Pa^\(1-components.count) s^-1"
        case .concentration(let concentration, _):
            value = concentration; standardUnits = "mol/L"; rateUnits = "(mol/L)^\(1-components.count) s^-1"
        }
        return try VivoConditionalEncounterPlan.prepare(request, components: components, sourceEndpoint: "reactants",
            temperatureK: 298.15, molecularity: molecularityOverride ?? components.count, standardState: standardState,
            standardStateValue: value, standardStateUnits: standardUnits,
            rateUnits: unitsOverride ?? (components.count == 1 ? "s^-1" : rateUnits), logRateConstant: logRate)
    }

    @Test func gasEncounterUsesBimolecularRateAndBothLogDerivatives() throws {
        let components = try components(), request = request(components)
        let prepared = try plan(request, components: components), result = try prepared.propagate(request)
        #expect(abs(prepared.hazard - 6) < 1e-13)
        #expect(abs(try #require(prepared.logHazard) - log(6)) < 1e-14)
        #expect(result.probabilityPoissonTailBound <= request.numerics.absoluteProbabilityTolerance)
        for observation in result.observations {
            let x = 6 * observation.timeSeconds, survival = exp(-x)
            #expect(abs(observation.survivalProbability - survival) < 2e-12)
            #expect(abs(observation.reactedProbability + expm1(-x)) < 2e-12)
            #expect(abs(observation.reactiveFluxPerSecond - 6 * survival) < 2e-11)
            #expect(abs(try #require(observation.hazardPerSecond) - 6) < 1e-12)
            #expect(observation.derivatives.map(\.parameterIdentifier) == ["natural-log-rate", "natural-log-reservoir-0"])
            for derivative in observation.derivatives {
                #expect(abs(derivative.survivalProbability + x * survival) < 2e-11)
                #expect(abs(derivative.reactedProbability - x * survival) < 2e-11)
                #expect(abs(derivative.reactiveFluxPerSecond - 6 * survival * (1 - x)) < 2e-10)
                #expect(abs(try #require(derivative.hazardPerSecond) - 6) < 2e-10)
            }
        }
    }

    @Test func concentrationProductsAndUnimolecularLimitRetainTheirOrders() throws {
        let three = try components(3)
        let request = request(three, values: [0.5, 0.25], unit: .molar)
        let prepared = try plan(request, components: three, logRate: log(8),
            standardState: .concentration(molPerLitre: 1, gasReferencePressurePa: 101325))
        #expect(abs(prepared.hazard - 1) < 1e-14)
        let result = try prepared.propagate(request)
        #expect(result.observations.last?.derivatives.map(\.parameterIdentifier) ==
            ["natural-log-rate", "natural-log-reservoir-0", "natural-log-reservoir-2"])
        #expect(abs(try #require(result.observations.last).survivalProbability - exp(-1)) < 2e-12)
        let one = try components(1), unimolecular = self.request(one, tagged: 0, values: [])
        let single = try plan(unimolecular, components: one)
        #expect(abs(single.hazard - 2) < 1e-14 && single.reservoirIndices.isEmpty)
        let observations = try single.propagate(unimolecular).observations
        #expect(observations.allSatisfy { $0.derivatives.count == 1 })
        #expect(abs(try #require(observations.last).survivalProbability - exp(-2)) < 2e-12)
    }

    @Test func zeroReservoirHasNoEventsAndDoesNotDefineALogReservoirDerivative() throws {
        let components = try components(), request = request(components, values: [0], sensitivities: false)
        let prepared = try plan(request, components: components)
        #expect(prepared.hazard == 0 && prepared.logHazard == nil)
        for observation in try prepared.propagate(request).observations {
            #expect(observation.survivalProbability == 1 && observation.reactedProbability == 0)
            #expect(observation.reactiveFluxPerSecond == 0 && observation.hazardPerSecond == 0)
            #expect(observation.derivatives.count == 1)
            #expect(observation.derivatives[0].survivalProbability == 0)
            #expect(observation.derivatives[0].hazardPerSecond == 0)
        }
        let invalid = self.request(components, values: [0], sensitivities: true)
        #expect(throws: (any Error).self) { try plan(invalid, components: components) }
    }

    @Test func unitsTemperatureAndMolecularityCannotBeRelabeled() throws {
        let components = try components()
        for invalid in [request(components, unit: .molar), request(components, temperature: 300),
                        request(components, endpoint: "products"), request(components, values: [-1]),
                        request(components, values: [.infinity]), request(components, values: [.nan])] {
            #expect(throws: (any Error).self) { try plan(invalid, components: components) }
        }
        #expect(throws: (any Error).self) { try plan(request(components), components: components, unitsOverride: "s^-1") }
        #expect(throws: (any Error).self) { try plan(request(components), components: components, molecularityOverride: 1) }
        #expect(throws: (any Error).self) {
            try plan(request(components, unit: .molar), components: components,
                standardState: .concentration(molPerLitre: 1, gasReferencePressurePa: 0))
        }
    }

    @Test func everyMappedComponentMustHaveOneExactRole() throws {
        let components = try components(3), good = request(components, values: [2, 3])
        func replace(_ tagged: VivoConditionalEncounterComponent, _ reservoirs: [VivoConditionalEncounterReservoir]) -> VivoConditionalEncounterRequest {
            .init(identifier: good.identifier, reactantEndpointIdentifier: good.reactantEndpointIdentifier,
                temperatureK: good.temperatureK, taggedComponent: tagged, reservoirs: reservoirs,
                observationTimesSeconds: good.observationTimesSeconds)
        }
        let changedPoint = VivoConditionalEncounterComponent(index: 0, atomIdentifiers: components[0].atomIdentifiers,
            qualifiedPointFingerprint: try VivoCanonicalJSON.fingerprint(Data("different qualified point".utf8)))
        let changedAtoms = VivoConditionalEncounterComponent(index: 0, atomIdentifiers: ["other atom"],
            qualifiedPointFingerprint: components[0].qualifiedPointFingerprint)
        let badIndex = VivoConditionalEncounterComponent(index: 99, atomIdentifiers: ["other atom"],
            qualifiedPointFingerprint: components[0].qualifiedPointFingerprint)
        let invalid = [
            replace(good.taggedComponent, Array(good.reservoirs.dropLast())),
            replace(good.taggedComponent, [good.reservoirs[0], good.reservoirs[0]]),
            replace(components[0], good.reservoirs),
            replace(badIndex, good.reservoirs),
            replace(good.taggedComponent, [.init(component: changedPoint, value: 2, unit: .pascal), good.reservoirs[1]]),
            replace(good.taggedComponent, [.init(component: changedAtoms, value: 2, unit: .pascal), good.reservoirs[1]])]
        for request in invalid { #expect(throws: (any Error).self) { try plan(request, components: components) } }
        let reordered = replace(good.taggedComponent, Array(good.reservoirs.reversed()))
        #expect(try plan(reordered, components: components).reservoirIndices == [0, 2])
    }

    @Test func logarithmicProductAvoidsIntermediateOverflowAndRejectsLostPositiveHazards() throws {
        let components = try components(3)
        let large = try plan(request(components, values: [1e200, 1e200]), components: components, logRate: log(1e-300))
        #expect(abs(large.hazard / 1e100 - 1) < 2e-12)
        let small = try plan(request(components, values: [1e-200, 1e-200]), components: components, logRate: log(1e300))
        #expect(abs(small.hazard / 1e-100 - 1) < 2e-12)
        let one = try self.components(1), request = self.request(one, tagged: 0, values: [])
        for logRate in [.nan, .infinity, log(Double.greatestFiniteMagnitude) + 1, log(Double.leastNonzeroMagnitude) - 1] {
            #expect(throws: (any Error).self) { try plan(request, components: one, logRate: logRate) }
        }
    }

    @Test func invalidSchedulesAndExistingProbabilityWorkLimitsReject() throws {
        let components = try components()
        for times: [Double] in [[], [-1], [0, 0], [1, 0], [.nan], [.infinity]] {
            #expect(throws: (any Error).self) { try plan(request(components, times: times), components: components) }
        }
        let enormous = request(components, times: [Double.greatestFiniteMagnitude])
        let prepared = try plan(enormous, components: components)
        #expect(throws: (any Error).self) { try prepared.propagate(enormous) }
        let tinyBudget = request(components, numerics: .init(maximumPrimitiveWork: 1))
        #expect(throws: (any Error).self) { try plan(tinyBudget, components: components).propagate(tinyBudget) }
    }

    /// Structurally valid mapped requests with deliberately absent path evidence.
    /// No electronic evaluation is performed: the authoritative validator must
    /// reject the absent four trials/twelve comparisons before any replay.
    private func unqualifiedSource() throws -> VivoTransitionStateTheoryResult {
        guard case .connectedReaction(let seed) = try VivoReactionQualificationWorkflow.template("h3-connected-rate").calculation else {
            throw VivoChemistryError.invalid("host reaction fixture template")
        }
        func point(_ request: VivoNuclearQualificationRequest) -> VivoNuclearQualifiedPoint {
            let positions = request.model.system.nuclei.map(\.positionBohr)
            let modes = VivoVibrationalAnalysis(positionsBohr: positions, massesDa: request.massesDa,
                cartesianHessian: VivoQMMatrix(0, 0), rigidModeCount: 0, signedFrequenciesCM: [],
                inertiaEigenvaluesDaBohr2: [], massWeightedModes: VivoQMMatrix(0, 0),
                maximumGradient: 0, gradientRMS: 0, rigidHessianResidual: 0)
            let thermo = VivoThermochemistryResult(configuration: request.thermochemistry, kind: request.kind,
                electronicEnergyHartree: 0, zeroPointEnergyHartree: 0, enthalpyCorrectionHartree: 0,
                gasEntropyHartreePerK: 0, standardStateShiftHartree: 0, gibbsEnergyHartree: 0,
                imaginaryFrequencyCM: nil, rotor: "unqualified host placeholder", modes: modes, model: "unqualified host placeholder")
            return .init(schema: VivoNuclearQualifiedPoint.schema, request: request, finalPositionsBohr: positions,
                thermochemistry: thermo, gradientStepAgreement: 0, hessianStepAgreement: 0, optimizationIterations: 0,
                electronicEvaluations: 0, qualification: "unqualified host placeholder")
        }
        let connection = VivoReactionConnectivityRequest(atomIdentifiers: seed.atomIdentifiers, saddle: point(seed.saddle),
            endpoints: seed.endpoints.map { endpoint in
                .init(identifier: endpoint.identifier, components: endpoint.components.map {
                    .init(atomIndices: $0.atomIndices, point: point($0.qualification))
                })
            }, configuration: seed.connectivity)
        let claimed = VivoReactionConnectivityResult(schema: VivoReactionConnectivityResult.schema, request: connection,
            trials: [], comparisons: [], descentElectronicEvaluations: 0, converged: true,
            interpretation: VivoReactionConnectivity.interpretation)
        return .init(schema: VivoTransitionStateTheoryResult.schema,
            request: .init(connectivity: claimed, reactantEndpointIdentifier: seed.reactantEndpointIdentifier,
                           transmission: seed.transmission),
            barrier: .init(electronicBarrierHartree: 0, zeroPointBarrierHartree: 0, activationGibbsHartree: 0,
                           standardStateContributionHartree: 0, reactantMolecularity: 2, model: "unqualified host placeholder"),
            temperatureK: 298.15, molecularity: 2, transmissionCoefficient: 1, standardStateValue: 101325,
            standardStateUnits: "Pa", rateUnits: "Pa^-1 s^-1", thermalFrequencyPerSecond: 1, activationExponent: 0,
            logRateConstant: log(2), rateConstant: 2, hasIndependentTransmissionEvidence: false,
            interpretation: VivoTransitionStateTheory.interpretation)
    }

    @Test func structuralBindingsAndZeroBathCannotPromoteADecodedQualificationClaim() throws {
        let source = try unqualifiedSource(), bindings = try VivoConditionalEncounter.componentBindings(in: source)
        #expect(bindings.map(\.atomIdentifiers) == [["H0", "H1"], ["H2"]])
        #expect(try VivoConditionalEncounter.componentBinding(index: 1, source: source) == bindings[1])
        #expect(throws: (any Error).self) { try VivoConditionalEncounter.componentBinding(index: 2, source: source) }
        for value in [0.0, 3.0] {
            let request = VivoConditionalEncounterRequest(identifier: "unqualified source must reject",
                reactantEndpointIdentifier: source.request.reactantEndpointIdentifier, temperatureK: source.temperatureK,
                taggedComponent: bindings[1], reservoirs: [.init(component: bindings[0], value: value, unit: .pascal)],
                observationTimesSeconds: [0, 1])
            try VivoConditionalEncounter.validateConditions(request, source: source)
            do {
                _ = try VivoConditionalEncounter.calculate(request, source: source)
                Issue.record("decoded convergence flag acquired full source authority")
            } catch VivoChemistryError.invalid(let message) {
                #expect(message == "reaction connectivity result identity or qualification contract")
            }
        }
    }

    @Test func workflowRejectsWrongReactionKindsAndUnqualifiedSourcesWithoutPublication() async throws {
        let source = try unqualifiedSource(), bindings = try VivoConditionalEncounter.componentBindings(in: source)
        let request = VivoConditionalEncounterRequest(identifier: "workflow source gate",
            reactantEndpointIdentifier: source.request.reactantEndpointIdentifier, temperatureK: source.temperatureK,
            taggedComponent: bindings[1], reservoirs: [.init(component: bindings[0], value: 0, unit: .pascal)],
            observationTimesSeconds: [0, 1])
        let requestData = try VivoCanonicalJSON.encode(request)
        let operation = VivoConditionalEncounterWorkflow.operation(implementationFingerprint:
            try VivoCanonicalJSON.fingerprint(Data("host implementation".utf8)))
        let wrong = try VivoCanonicalJSON.encode(VivoReactionCalculationResult.connectivity(result: source.request.connectivity))
        let unqualified = try VivoCanonicalJSON.encode(VivoReactionCalculationResult.connectedReaction(result: source))
        for data in [wrong, unqualified] {
            do {
                _ = try await operation.execute(.object([:]), ["request": requestData, "reaction": data], .init())
                Issue.record("invalid source yielded workflow outputs")
            } catch VivoChemistryError.invalid(let message) {
                #expect(message == "conditional encounter needs a complete molecular TST reaction result" ||
                        message == "reaction connectivity result identity or qualification contract")
            }
        }
        do {
            _ = try await operation.execute(.object([:]), ["request": requestData, "reaction": unqualified], .init(maximumBytes: 1))
            Issue.record("aggregate input budget did not reject")
        } catch VivoChemistryError.invalid(let message) {
            #expect(message == "conditional encounter input slots or aggregate byte budget")
        }
        #expect(throws: (any Error).self) {
            try operation.validateOutputs(.object([:]), ["request": requestData, "reaction": unqualified], [:], .init())
        }
    }
}
