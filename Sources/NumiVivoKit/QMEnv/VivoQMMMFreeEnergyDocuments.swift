import Foundation

public struct VivoQMMMFreeEnergyAnalysisRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-analysis/v1"
    public var schema:String
    public var coordinate:VivoQMMMReactionCoordinate
    public var temperatureK:Double
    public var traces:[VivoQMMMUmbrellaTrace]
    public var configuration:VivoQMMMFreeEnergyAnalysisConfiguration
    public init(coordinate:VivoQMMMReactionCoordinate,temperatureK:Double,traces:[VivoQMMMUmbrellaTrace],configuration:VivoQMMMFreeEnergyAnalysisConfiguration) {
        schema=Self.schema;self.coordinate=coordinate;self.temperatureK=temperatureK;self.traces=traces;self.configuration=configuration
    }
    public func calculate() throws -> VivoQMMMActivationFreeEnergyResult {
        guard schema==Self.schema else { throw VivoChemistryError.invalid("QM/MM free-energy analysis schema") }
        return try VivoQMMMFreeEnergy.analyze(coordinate:coordinate,temperatureK:temperatureK,traces:traces,configuration:configuration)
    }
}
