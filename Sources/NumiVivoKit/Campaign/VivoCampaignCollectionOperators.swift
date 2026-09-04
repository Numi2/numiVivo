import Foundation

/// Exact comparison used by the deterministic campaign bin-packer when the
/// policy requires a single numerical authority per batch.
func == (
    lhs: Set<VivoAdaptiveFidelityMode>,
    rhs: [VivoAdaptiveFidelityMode]
) -> Bool {
    lhs == Set(rhs)
}

func == (
    lhs: [VivoAdaptiveFidelityMode],
    rhs: Set<VivoAdaptiveFidelityMode>
) -> Bool {
    Set(lhs) == rhs
}
