import Foundation

/// Principal polarizabilities use the classical geometry unit nm^3. The
/// electronic operator converts explicitly to Bohr^3 at its boundary.
public struct VivoInducedDipoleSite: Codable, Sendable, Equatable {
    public let particleIndex: UInt32
    public let principalPolarizabilitiesNM3: VivoVector3D
    /// z points from the polarizable center to zParticle; x is the orthogonal
    /// projection toward xParticle. An anisotropic site requires this frame.
    public let zParticle: UInt32?
    public let xParticle: UInt32?
    public init(particleIndex: UInt32,principalPolarizabilitiesNM3: VivoVector3D,zParticle: UInt32? = nil,xParticle: UInt32? = nil) {
        self.particleIndex=particleIndex;self.principalPolarizabilitiesNM3=principalPolarizabilitiesNM3
        self.zParticle=zParticle;self.xParticle=xParticle
    }
}
public struct VivoPolarizationPair: Codable, Sendable, Equatable {
    public let firstParticle: UInt32
    public let secondParticle: UInt32
    /// Exponential Thole v=screening*r^3 for the declared primary pair. nil
    /// leaves the radial kernel undamped while retaining explicit pair scales.
    public let screeningPerNM3: Double?
    public let permanentFieldScale: Double
    public let mutualInductionScale: Double
    public init(firstParticle: UInt32,secondParticle: UInt32,screeningPerNM3: Double? = nil,
                permanentFieldScale: Double = 1,mutualInductionScale: Double = 1) {
        self.firstParticle=firstParticle;self.secondParticle=secondParticle;self.screeningPerNM3=screeningPerNM3
        self.permanentFieldScale=permanentFieldScale;self.mutualInductionScale=mutualInductionScale
    }
}
public struct VivoInducedDipoleConfiguration: Codable, Sendable, Equatable {
    public static let model = "variational-induced-dipoles-exp-Thole-primary-pairs-v1"
    public let model: String
    public var sites: [VivoInducedDipoleSite]
    public var pairs: [VivoPolarizationPair]
    public var residualToleranceHartreePerEBohr: Double
    public var minimumResponsePivot: Double
    public var maximumPolarizableSites: Int
    public var parameterProvenance: String
    public init(sites: [VivoInducedDipoleSite],pairs: [VivoPolarizationPair] = [],
                residualToleranceHartreePerEBohr: Double = 1e-10,minimumResponsePivot: Double = 1e-10,
                maximumPolarizableSites: Int = 64,parameterProvenance: String) {
        model=Self.model;self.sites=sites;self.pairs=pairs
        self.residualToleranceHartreePerEBohr=residualToleranceHartreePerEBohr;self.minimumResponsePivot=minimumResponsePivot
        self.maximumPolarizableSites=maximumPolarizableSites;self.parameterProvenance=parameterProvenance
    }
    public func validate(particleCount: Int? = nil) throws {
        guard model==Self.model,!sites.isEmpty,Set(sites.map(\.particleIndex)).count==sites.count,
              (1...256).contains(maximumPolarizableSites),sites.count<=maximumPolarizableSites,!parameterProvenance.isEmpty,
              residualToleranceHartreePerEBohr.isFinite,residualToleranceHartreePerEBohr>0,
              minimumResponsePivot.isFinite,minimumResponsePivot>0 else { throw VivoChemistryError.invalid("induced-dipole model, provenance or numerical limits") }
        var ids = Set<UInt64>()
        func index(_ i: UInt32) -> Bool { particleCount == nil || Int(i)<particleCount! }
        for site in sites {
            let alpha=site.principalPolarizabilitiesNM3
            guard index(site.particleIndex),alpha.isFinite,alpha.x>0,alpha.y>0,alpha.z>0,
                  (site.zParticle == nil) == (site.xParticle == nil) else { throw VivoChemistryError.invalid("polarizability tensor or frame ownership") }
            if let z=site.zParticle,let x=site.xParticle {
                guard index(z),index(x),Set([site.particleIndex,z,x]).count==3 else { throw VivoChemistryError.invalid("polarization frame identities") }
            } else {
                guard alpha.x==alpha.y,alpha.x==alpha.z else { throw VivoChemistryError.invalid("anisotropic polarizability requires a differentiable local frame") }
            }
        }
        for pair in pairs {
            let a=min(pair.firstParticle,pair.secondParticle),b=max(pair.firstParticle,pair.secondParticle)
            guard index(a),index(b),a != b,ids.insert(UInt64(a)<<32|UInt64(b)).inserted,
                  pair.permanentFieldScale.isFinite,pair.permanentFieldScale>=0,
                  pair.mutualInductionScale.isFinite,pair.mutualInductionScale>=0,
                  pair.screeningPerNM3 == nil || (pair.screeningPerNM3!.isFinite && pair.screeningPerNM3!>0) else {
                throw VivoChemistryError.invalid("polarization pair damping or scaling")
            }
        }
    }
}
public struct VivoInducedDipoleMoment: Codable, Sendable, Equatable {
    public let particleIndex: UInt32
    public let dipoleEBohr: VivoVector3D
}
public struct VivoInducedDipoleResult: Codable, Sendable, Equatable {
    public let model: String
    public let moments: [VivoInducedDipoleMoment]
    public let maximumResponseResidual: Double
    public let minimumResponsePivot: Double
    public let selfEnergyHartree: Double
    public let dampingCorrectionEnergyHartree: Double
    public let physicalChargeChannelCorrectionHartree: Double
    /// Physical MM charge channel in the electronic source point-charge order.
    /// It may differ from the boundary-transformed embedding charges.
    public let physicalPermanentChargesE: [Double]
}
