import Foundation

#if canImport(CoreML)
@preconcurrency import CoreML

/// Apple-native Core ML backend for compiled surrogate blocks. Model outputs are
/// still passed through `VivoSurrogateAuthorityGate`; Core ML never owns or
/// publishes authoritative biological state.
public actor VivoCoreMLSurrogateBackend: VivoSurrogateBackend {
    public nonisolated let backendID = "numivivo.coreml.v1"
    public nonisolated let modelFingerprint: String
    public nonisolated let inputWidth: Int
    public nonisolated let outputWidth: Int

    private let model: MLModel
    private let inputFeatureName: String
    private let outputFeatureName: String
    private let uncertaintyFeatureName: String?
    private let fixedUncertainty: [Float]?

    public init(
        modelURL: URL,
        modelFingerprint: String,
        inputWidth: Int,
        outputWidth: Int,
        inputFeatureName: String = "input",
        outputFeatureName: String = "output",
        uncertaintyFeatureName: String? = nil,
        fixedUncertainty: [Float]? = nil,
        computeUnits: MLComputeUnits = .all
    ) throws {
        guard !modelFingerprint.isEmpty,
              inputWidth > 0,
              outputWidth > 0,
              !inputFeatureName.isEmpty,
              !outputFeatureName.isEmpty else {
            throw VivoSurrogateError.backendContractMismatch
        }
        if let fixedUncertainty {
            guard fixedUncertainty.count == outputWidth,
                  fixedUncertainty.allSatisfy({ $0.isFinite && $0 >= 0 }) else {
                throw VivoSurrogateError.backendContractMismatch
            }
        }

        let configuration = MLModelConfiguration()
        configuration.computeUnits = computeUnits
        do {
            self.model = try MLModel(contentsOf: modelURL, configuration: configuration)
        } catch {
            throw VivoSurrogateError.backendUnavailable(error.localizedDescription)
        }
        self.modelFingerprint = modelFingerprint
        self.inputWidth = inputWidth
        self.outputWidth = outputWidth
        self.inputFeatureName = inputFeatureName
        self.outputFeatureName = outputFeatureName
        self.uncertaintyFeatureName = uncertaintyFeatureName
        self.fixedUncertainty = fixedUncertainty
    }

    public func predict(_ batch: VivoSurrogateBatch) async throws -> VivoSurrogatePrediction {
        guard batch.inputWidth == inputWidth else {
            throw VivoSurrogateError.invalidBatch
        }
        let input: MLMultiArray
        do {
            input = try MLMultiArray(
                shape: [NSNumber(value: batch.batchSize), NSNumber(value: inputWidth)],
                dataType: .float32
            )
        } catch {
            throw VivoSurrogateError.modelFailure("input allocation: \(error.localizedDescription)")
        }
        for index in batch.values.indices {
            input[index] = NSNumber(value: batch.values[index])
        }

        let provider: MLDictionaryFeatureProvider
        do {
            provider = try MLDictionaryFeatureProvider(dictionary: [
                inputFeatureName: MLFeatureValue(multiArray: input)
            ])
        } catch {
            throw VivoSurrogateError.modelFailure("input provider: \(error.localizedDescription)")
        }

        let result: MLFeatureProvider
        do {
            result = try model.prediction(from: provider)
        } catch {
            throw VivoSurrogateError.modelFailure(error.localizedDescription)
        }
        guard let output = result.featureValue(for: outputFeatureName)?.multiArrayValue else {
            throw VivoSurrogateError.modelFailure("missing output feature \(outputFeatureName)")
        }
        let expected = batch.batchSize.multipliedReportingOverflow(by: outputWidth)
        guard !expected.overflow, output.count == expected.partialValue else {
            throw VivoSurrogateError.modelFailure(
                "output feature has \(output.count) elements; expected \(expected.partialValue)"
            )
        }
        var values = [Float](repeating: 0, count: output.count)
        for index in values.indices {
            values[index] = output[index].floatValue
        }

        let uncertainty: [Float]
        if let uncertaintyFeatureName {
            guard let array = result.featureValue(for: uncertaintyFeatureName)?.multiArrayValue else {
                throw VivoSurrogateError.modelFailure("missing uncertainty feature \(uncertaintyFeatureName)")
            }
            guard array.count == batch.batchSize || array.count == expected.partialValue else {
                throw VivoSurrogateError.modelFailure("uncertainty output has an invalid shape")
            }
            uncertainty = (0..<array.count).map { array[$0].floatValue }
        } else if let fixedUncertainty {
            var expanded = [Float](repeating: 0, count: expected.partialValue)
            for row in 0..<batch.batchSize {
                for column in 0..<outputWidth {
                    expanded[row * outputWidth + column] = fixedUncertainty[column]
                }
            }
            uncertainty = expanded
        } else {
            uncertainty = [Float](repeating: 0, count: batch.batchSize)
        }

        return try .init(
            batchSize: batch.batchSize,
            outputWidth: outputWidth,
            values: values,
            normalizedUncertainty: uncertainty,
            backend: backendID,
            modelFingerprint: modelFingerprint
        )
    }
}
#endif
