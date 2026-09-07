import Foundation

/// Compatibility surface for platform metric prose. The replicated-rate result
/// intentionally stores numerical evidence rather than duplicated prose; platform
/// adapters may still expose a stable interpretation without changing its schema.
public extension VivoQMMMReplicatedFreeEnergyRateResult {
    var interpretation: String {
        "Replicated conditional QM/MM rate evidence with explicit between-replica dispersion and unresolved within-replica uncertainty preserved as nil. This result does not include missing model, mechanism, state-population or experimental uncertainty."
    }
}
