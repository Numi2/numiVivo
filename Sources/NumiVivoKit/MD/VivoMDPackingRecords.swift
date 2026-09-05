// The packer and its serialized GPU records deliberately have different owners.
// Resolve short names within the packer's scope to the existing record types;
// never introduce parallel structs whose memory layout could drift from Metal.
extension VivoMDSystemPacker {
    typealias Bond = VivoMDPackedSystem.Bond
    typealias Angle = VivoMDPackedSystem.Angle
    typealias Torsion = VivoMDPackedSystem.Torsion
    typealias Constraint = VivoMDPackedSystem.Constraint
    typealias PairException = VivoMDPackedSystem.PairException
}
