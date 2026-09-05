import Foundation

public enum VivoUnimolecularEyring {
    /// Scalar unimolecular Eyring relation. The argument MUST be a Gibbs
    /// activation energy, not a raw electronic-energy difference. Use
    /// VivoReactionEnergies for provenance/reference/stationary-point checks.
    public static func logRatePerSecond(activationGibbsHartree: Double, temperatureK: Double,
                                        transmissionCoefficient: Double = 1) throws -> Double {
        guard activationGibbsHartree.isFinite,temperatureK.isFinite,temperatureK>0,
              transmissionCoefficient.isFinite,transmissionCoefficient>0,transmissionCoefficient<=1 else {
            throw VivoChemistryError.invalid("Eyring Gibbs barrier, temperature or transmission coefficient")
        }
        let exponent=activationGibbsHartree*VivoAtomicUnits.hartreeInKJPerMol*1000/(VivoAtomicUnits.gasConstantJPerMolK*temperatureK)
        let result=log(transmissionCoefficient)+log(VivoAtomicUnits.boltzmannJPerK*temperatureK/VivoAtomicUnits.planckJS)-exponent
        guard result.isFinite else { throw VivoChemistryError.invalid("Eyring log-rate overflow") };return result
    }
}
