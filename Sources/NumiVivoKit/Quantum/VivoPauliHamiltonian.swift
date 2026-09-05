import Foundation

public struct VivoComplex: Codable, Sendable, Equatable, Hashable {
    public var real: Double
    public var imaginary: Double
    public init(_ real: Double, _ imaginary: Double = 0) { self.real = real; self.imaginary = imaginary }
    public static let zero = VivoComplex(0)
    public var magnitudeSquared: Double { real*real + imaginary*imaginary }
    public var conjugate: Self { .init(real,-imaginary) }
    public static func +(lhs: Self, rhs: Self) -> Self { .init(lhs.real+rhs.real,lhs.imaginary+rhs.imaginary) }
    public static func -(lhs: Self, rhs: Self) -> Self { .init(lhs.real-rhs.real,lhs.imaginary-rhs.imaginary) }
    public static func *(lhs: Self, rhs: Self) -> Self { .init(lhs.real*rhs.real-lhs.imaginary*rhs.imaginary,lhs.real*rhs.imaginary+lhs.imaginary*rhs.real) }
    public static func *(lhs: Self, rhs: Double) -> Self { .init(lhs.real*rhs,lhs.imaginary*rhs) }
}

public struct VivoPauliKey: Codable, Sendable, Equatable, Hashable, Comparable {
    public let xMask: UInt64
    public let zMask: UInt64
    public init(xMask: UInt64, zMask: UInt64) { self.xMask=xMask; self.zMask=zMask }
    public static func <(lhs: Self, rhs: Self) -> Bool { lhs.xMask == rhs.xMask ? lhs.zMask < rhs.zMask : lhs.xMask < rhs.xMask }
}
public struct VivoPauliTerm: Codable, Sendable, Equatable {
    public let key: VivoPauliKey
    public let coefficient: VivoComplex
}
public struct VivoPauliHamiltonian: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/pauli-hamiltonian/v1"
    public let schema: String
    public let qubitCount: Int
    public let terms: [VivoPauliTerm]
    public let sourceEnergyReference: String
    public let provenance: [String:String]
    public init(qubitCount: Int, terms: [VivoPauliTerm], sourceEnergyReference: String,
                provenance: [String:String] = [:]) {
        schema=Self.schema; self.qubitCount=qubitCount; self.terms=terms; self.sourceEnergyReference=sourceEnergyReference; self.provenance=provenance
    }
    public func validate(tolerance: Double = 1e-12) throws {
        guard schema==Self.schema,(1...62).contains(qubitCount),!sourceEnergyReference.isEmpty,
              tolerance.isFinite,tolerance>0 else { throw VivoChemistryError.invalid("Pauli Hamiltonian metadata") }
        let validMask=(UInt64(1)<<qubitCount)-1
        var seen=Set<VivoPauliKey>()
        for term in terms {
            guard term.key.xMask & ~validMask == 0,term.key.zMask & ~validMask == 0,
                  term.coefficient.real.isFinite,term.coefficient.imaginary.isFinite,
                  seen.insert(term.key).inserted else { throw VivoChemistryError.invalid("Pauli term mask/coefficient/duplication") }
            guard abs(term.coefficient.imaginary)<=tolerance else { throw VivoChemistryError.invalid("Hermitian electronic Hamiltonian produced an imaginary Pauli coefficient") }
        }
    }
}

private enum VivoLocalPauli: Int { case i=0,x=1,y=2,z=3 }
private func vivoLocalPauli(_ x: Bool,_ z: Bool) -> VivoLocalPauli {
    if x { return z ? .y : .x }; return z ? .z : .i
}
private func vivoPauliBits(_ p: VivoLocalPauli) -> (Bool,Bool) {
    switch p { case .i:return(false,false);case .x:return(true,false);case .y:return(true,true);case .z:return(false,true) }
}
private func vivoMultiplyLocal(_ a: VivoLocalPauli,_ b: VivoLocalPauli) -> (VivoLocalPauli,VivoComplex) {
    if a == .i { return (b,.init(1)) }; if b == .i { return (a,.init(1)) }; if a == b { return (.i,.init(1)) }
    switch (a,b) {
    case (.x,.y): return (.z,.init(0,1)); case (.y,.x): return (.z,.init(0,-1))
    case (.y,.z): return (.x,.init(0,1)); case (.z,.y): return (.x,.init(0,-1))
    case (.z,.x): return (.y,.init(0,1)); case (.x,.z): return (.y,.init(0,-1))
    default: fatalError("exhaustive Pauli multiplication")
    }
}
private func vivoMultiplyPauli(_ a: VivoPauliKey,_ b: VivoPauliKey,qubits: Int) -> (VivoPauliKey,VivoComplex) {
    var x:UInt64=0,z:UInt64=0,phase=VivoComplex(1)
    for q in 0..<qubits {
        let bit=UInt64(1)<<q
        let pa=vivoLocalPauli(a.xMask&bit != 0,a.zMask&bit != 0)
        let pb=vivoLocalPauli(b.xMask&bit != 0,b.zMask&bit != 0)
        let (pc,p)=vivoMultiplyLocal(pa,pb); phase=phase*p
        let bits=vivoPauliBits(pc); if bits.0 { x|=bit }; if bits.1 { z|=bit }
    }
    return (.init(xMask:x,zMask:z),phase)
}

private struct VivoJWComponent { let key: VivoPauliKey; let coefficient: VivoComplex }
private enum VivoJWOperator { case creation(Int), annihilation(Int) }
private func vivoJW(_ op: VivoJWOperator) -> [VivoJWComponent] {
    let mode:Int,creation:Bool
    switch op { case .creation(let m):mode=m;creation=true;case .annihilation(let m):mode=m;creation=false }
    let parityMask = mode == 0 ? UInt64(0) : (UInt64(1)<<mode)-1
    let bit=UInt64(1)<<mode
    let x=VivoPauliKey(xMask:bit,zMask:parityMask)
    let y=VivoPauliKey(xMask:bit,zMask:parityMask|bit)
    return [.init(key:x,coefficient:.init(0.5)),.init(key:y,coefficient:.init(0,creation ? -0.5 : 0.5))]
}

public enum VivoJordanWigner {
    private static func expand(_ operators: [VivoJWOperator], coefficient: Double, qubits: Int,
                               into accumulator: inout [VivoPauliKey:VivoComplex]) {
        var components=[VivoJWComponent(key:.init(xMask:0,zMask:0),coefficient:.init(coefficient))]
        for op in operators {
            let mapped=vivoJW(op)
            var next:[VivoJWComponent]=[]; next.reserveCapacity(components.count*2)
            for left in components { for right in mapped {
                let product=vivoMultiplyPauli(left.key,right.key,qubits:qubits)
                next.append(.init(key:product.0,coefficient:left.coefficient*right.coefficient*product.1))
            } }
            components=next
        }
        for term in components { accumulator[term.key]=(accumulator[term.key] ?? .zero)+term.coefficient }
    }
    public static func transform(_ h: VivoEmbeddedHamiltonian, coefficientTolerance: Double = 1e-12,
                                 budget: VivoChemistryBudget = .init()) throws -> VivoPauliHamiltonian {
        try h.validate(budget:budget)
        guard coefficientTolerance.isFinite,coefficientTolerance>0 else { throw VivoChemistryError.invalid("JW coefficient tolerance") }
        let n=h.orbitalCount,qubits=2*n
        guard qubits<=62 else { throw VivoChemistryError.resourceLimit("Jordan-Wigner mask width") }
        var terms:[VivoPauliKey:VivoComplex]=[.init(xMask:0,zMask:0):.init(h.constantEnergyHartree)]
        for p in 0..<n { for q in 0..<n where abs(h.oneElectron[p,q])>coefficientTolerance {
            for spin in 0..<2 {
                expand([.creation(2*p+spin),.annihilation(2*q+spin)],coefficient:h.oneElectron[p,q],qubits:qubits,into:&terms)
            }
        } }
        for p in 0..<n { for q in 0..<n { for r in 0..<n { for s in 0..<n {
            let value=0.5*h.eri(p,q,r,s); if abs(value)<=coefficientTolerance { continue }
            for sigma in 0..<2 { for tau in 0..<2 {
                expand([.creation(2*p+sigma),.creation(2*r+tau),.annihilation(2*s+tau),.annihilation(2*q+sigma)],
                       coefficient:value,qubits:qubits,into:&terms)
            } }
        } } } }
        let cleaned=terms.filter{$0.value.magnitudeSquared>coefficientTolerance*coefficientTolerance}.map{VivoPauliTerm(key:$0.key,coefficient:$0.value)}.sorted{$0.key<$1.key}
        let result=VivoPauliHamiltonian(qubitCount:qubits,terms:cleaned,sourceEnergyReference:h.energyReference,
            provenance:h.provenance.merging(["mapping":"Jordan-Wigner","spinOrder":"interleaved-alpha-beta","coefficientTolerance":String(coefficientTolerance)]){$1})
        try result.validate(tolerance:max(1e-10,10*coefficientTolerance)); return result
    }
}
