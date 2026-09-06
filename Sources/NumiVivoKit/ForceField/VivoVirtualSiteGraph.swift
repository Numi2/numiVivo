import Foundation

public enum VivoSiteConstructionRule: Codable, Sendable, Equatable {
    case linear(weights: [Double])
    /// r0 + a*(r1-r0) + b*(r2-r0) + c*((r1-r0) x (r2-r0)).
    case outOfPlane(weight12: Double, weight13: Double, crossWeightPerNM: Double)
    case localCoordinates(originWeights: [Double], xWeights: [Double], yWeights: [Double], positionNM: VivoVector3D)
}
public struct VivoDependentSite: Codable, Sendable, Equatable {
    public let siteParticle: UInt32
    public let parentParticles: [UInt32]
    public let rule: VivoSiteConstructionRule
    public let provenance: String
    public init(siteParticle: UInt32, parentParticles: [UInt32], rule: VivoSiteConstructionRule, provenance: String) {
        self.siteParticle = siteParticle; self.parentParticles = parentParticles; self.rule = rule; self.provenance = provenance
    }
    public init(_ legacy: VivoLinearVirtualSite) {
        self.init(siteParticle: legacy.siteParticle, parentParticles: legacy.parentParticles,
                  rule: .linear(weights: legacy.weights), provenance: legacy.provenance ?? "legacy-linear-site-v1")
    }
}
public struct VivoSiteJacobian: Codable, Sendable, Equatable {
    /// Cartesian columns dr_site / d(parent.x,parent.y,parent.z).
    public let x: VivoVector3D
    public let y: VivoVector3D
    public let z: VivoVector3D
    public func apply(_ v: VivoVector3D) -> VivoVector3D { x*v.x+y*v.y+z*v.z }
    public func transposeApply(_ force: VivoVector3D) -> VivoVector3D { .init(x.dot(force),y.dot(force),z.dot(force)) }
}
public struct VivoDependentSiteState: Sendable {
    public let positionsNM: [VivoVector3D]
    /// One Jacobian per direct parent in graph site/parent order.
    public let jacobians: [[VivoSiteJacobian]]
}

/// One topologically ordered construction and reverse-force authority. All sites
/// retain their original particle IDs; parent dependencies never become atoms.
public struct VivoVirtualSiteGraph: Codable, Sendable, Equatable {
    public let particleCount: Int
    public let physicalParticles: [UInt32]
    public let sites: [VivoDependentSite]
    public let depths: [Int]
    public let physicalAncestors: [[UInt32]]
    public var maximumDepth: Int { depths.max() ?? -1 }

    public init(particleCount: Int, physicalParticles: [UInt32], definitions: [VivoDependentSite]) throws {
        guard particleCount > 0, particleCount <= Int(UInt32.max),
              Set(physicalParticles).count == physicalParticles.count,
              physicalParticles.allSatisfy({ Int($0) < particleCount }),
              Set(definitions.map(\.siteParticle)).count == definitions.count,
              Set(physicalParticles).isDisjoint(with: Set(definitions.map(\.siteParticle))),
              physicalParticles.count + definitions.count == particleCount else {
            throw VivoChemistryError.invalid("dependent-site graph particle ownership")
        }
        for site in definitions {
            guard Int(site.siteParticle) < particleCount, (1...4).contains(site.parentParticles.count),
                  !site.provenance.isEmpty, Set(site.parentParticles).count == site.parentParticles.count,
                  site.parentParticles.allSatisfy({ Int($0) < particleCount && $0 != site.siteParticle }) else {
                throw VivoChemistryError.invalid("dependent-site indices, parents or provenance")
            }
            func weights(_ w: [Double], sum: Double) -> Bool {
                w.count == site.parentParticles.count && w.allSatisfy(\.isFinite) && abs(w.reduce(0,+)-sum) <= 1e-9
            }
            switch site.rule {
            case .linear(let w):
                guard weights(w,sum: 1) else { throw VivoChemistryError.invalid("linear-site weights must sum to one") }
            case .outOfPlane(let a,let b,let c):
                guard site.parentParticles.count == 3, [a,b,c].allSatisfy(\.isFinite) else {
                    throw VivoChemistryError.invalid("out-of-plane site parameters")
                }
            case .localCoordinates(let o,let x,let y,let p):
                guard weights(o,sum: 1),weights(x,sum: 0),weights(y,sum: 0),p.isFinite,
                      x.contains(where: { $0 != 0 }),y.contains(where: { $0 != 0 }) else {
                    throw VivoChemistryError.invalid("local-frame site affine weights or local position")
                }
            }
        }
        var resolved = Dictionary(uniqueKeysWithValues: physicalParticles.map { ($0,(-1,[$0])) })
        var pending = definitions.sorted { $0.siteParticle < $1.siteParticle }, ordered: [VivoDependentSite] = []
        var levels: [Int] = [], ancestors: [[UInt32]] = []
        while !pending.isEmpty {
            let ready = pending.filter { site in site.parentParticles.allSatisfy { resolved[$0] != nil } }
            guard !ready.isEmpty else { throw VivoChemistryError.invalid("cyclic or unresolved dependent-site graph") }
            let ids = Set(ready.map(\.siteParticle))
            for site in ready {
                let parents = site.parentParticles.map { resolved[$0]! }
                let depth = (parents.map(\.0).max() ?? -1)+1
                let leaves = Set(parents.flatMap(\.1)).sorted()
                resolved[site.siteParticle] = (depth,leaves)
                ordered.append(site); levels.append(depth); ancestors.append(leaves)
            }
            pending.removeAll { ids.contains($0.siteParticle) }
        }
        self.particleCount = particleCount; self.physicalParticles = physicalParticles.sorted()
        sites = ordered; depths = levels; physicalAncestors = ancestors
    }
    public func validate() throws {
        guard try Self(particleCount: particleCount,physicalParticles: physicalParticles,definitions: sites) == self else {
            throw VivoChemistryError.invalid("decoded dependent-site graph does not reconstruct")
        }
    }
    /// cell=nil consumes finite, coherent molecular coordinates. With a cell,
    /// local parent displacements are made coherent before evaluating the rule.
    public func construct(positionsNM input: [VivoVector3D], periodicCell: VivoPeriodicCell? = nil) throws -> VivoDependentSiteState {
        try validate()
        guard input.count == particleCount,input.allSatisfy(\.isFinite),periodicCell?.isValid != false else {
            throw VivoChemistryError.invalid("dependent-site coordinate shape or cell")
        }
        var positions = input, jacobians: [[VivoSiteJacobian]] = []
        for site in sites {
            let base = positions[Int(site.parentParticles[0])]
            let parents = try site.parentParticles.map { index -> VivoVector3D in
                let p = positions[Int(index)]
                return try periodicCell.map { base + (try $0.minimumImage(p-base)) } ?? p
            }
            let state = try Self.construct(site.rule,parents: parents)
            positions[Int(site.siteParticle)] = state.position; jacobians.append(state.jacobians)
        }
        return .init(positionsNM: positions,jacobians: jacobians)
    }
    public func velocities(_ physicalVelocities: [VivoVector3D], state: VivoDependentSiteState) throws -> [VivoVector3D] {
        try validateState(state,values: physicalVelocities)
        var velocities = physicalVelocities
        for (i,site) in sites.enumerated() {
            velocities[Int(site.siteParticle)] = zip(site.parentParticles,state.jacobians[i]).reduce(.zero) {
                $0+$1.1.apply(velocities[Int($1.0)])
            }
        }
        guard velocities.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("dependent-site velocity overflow") }
        return velocities
    }
    /// Input forces must be RAW center forces. Calling this twice is not valid.
    public func redistribute(rawForces: [VivoVector3D], state: VivoDependentSiteState) throws -> [VivoVector3D] {
        try validateState(state,values: rawForces)
        var forces = rawForces
        for i in sites.indices.reversed() {
            let site = sites[i], f = forces[Int(site.siteParticle)]
            for (parent,jacobian) in zip(site.parentParticles,state.jacobians[i]) {
                forces[Int(parent)] = forces[Int(parent)]+jacobian.transposeApply(f)
            }
            forces[Int(site.siteParticle)] = .zero
        }
        guard forces.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("dependent-site force overflow") }
        return forces
    }
    /// Extra dE/d(strain) from nonlinear/nested construction, relative to
    /// independently affine-moving site centers. Fixed lattice image coefficients
    /// are differentiated with the cell; wrapped parent coordinates are not used
    /// as substitutes for their coherent local images.
    public func affineStrainCorrection(rawForces: [VivoVector3D],state: VivoDependentSiteState,
                                       periodicCell: VivoPeriodicCell? = nil,
                                       positionUnitScale: Double = 1) throws -> VivoQMMatrix {
        try validateState(state,values:rawForces)
        guard positionUnitScale.isFinite,positionUnitScale>0,periodicCell?.isValid != false else {
            throw VivoChemistryError.invalid("site strain units or cell")
        }
        var result=VivoQMMatrix(3,3)
        func affine(_ v: VivoVector3D,_ a: Int,_ b: Int) -> VivoVector3D {
            let x=[v.x,v.y,v.z][b]
            return a==0 ? .init(x,0,0):(a==1 ? .init(0,x,0):.init(0,0,x))
        }
        for a in 0..<3 { for b in 0..<3 {
            var tangent=state.positionsNM.map { affine($0,a,b) }
            for (i,site) in sites.enumerated() {
                let origin=state.positionsNM[Int(site.parentParticles[0])]
                var velocity=VivoVector3D.zero
                for (parent,jacobian) in zip(site.parentParticles,state.jacobians[i]) {
                    let p=state.positionsNM[Int(parent)]
                    let coherent=try periodicCell.map { origin+(try $0.minimumImage(p-origin)) } ?? p
                    velocity=velocity+jacobian.apply(tangent[Int(parent)]+affine(coherent-p,a,b))
                }
                tangent[Int(site.siteParticle)]=velocity
            }
            for site in sites {
                let i=Int(site.siteParticle),nonAffine=tangent[i]-affine(state.positionsNM[i],a,b)
                result[a,b]-=rawForces[i].dot(nonAffine)*positionUnitScale
            }
        } }
        guard result.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("site strain derivative overflow") }
        return result
    }
    private func validateState(_ state: VivoDependentSiteState, values: [VivoVector3D]) throws {
        guard state.positionsNM.count == particleCount, values.count == particleCount, values.allSatisfy(\.isFinite),
              state.jacobians.count == sites.count,
              zip(sites,state.jacobians).allSatisfy({ $0.parentParticles.count == $1.count }) else {
            throw VivoChemistryError.invalid("dependent-site derivative state shape")
        }
    }
    static func construct(_ rule: VivoSiteConstructionRule, parents p: [VivoVector3D]) throws -> (position: VivoVector3D,jacobians: [VivoSiteJacobian]) {
        let axes: [VivoVector3D] = [.init(1,0,0),.init(0,1,0),.init(0,0,1)]
        var position = VivoVector3D.zero, columns = [[VivoVector3D]](repeating: Array(repeating: .zero,count: 3),count: p.count)
        switch rule {
        case .linear(let weights):
            for i in p.indices {
                position = position+p[i]*weights[i]
                for a in 0..<3 { columns[i][a] = axes[a]*weights[i] }
            }
        case .outOfPlane(let a,let b,let c):
            let u = p[1]-p[0],v = p[2]-p[0]
            position = p[0]+u*a+v*b+u.cross(v)*c
            for i in p.indices { for axis in 0..<3 {
                let d0 = i == 0 ? axes[axis] : .zero, d1 = i == 1 ? axes[axis] : .zero, d2 = i == 2 ? axes[axis] : .zero
                let du = d1-d0,dv = d2-d0
                columns[i][axis] = d0+du*a+dv*b+(du.cross(v)+u.cross(dv))*c
            } }
        case .localCoordinates(let ow,let xw,let yw,let local):
            func weighted(_ weights: [Double]) -> VivoVector3D { zip(p,weights).reduce(.zero) { $0+$1.0*$1.1 } }
            let origin = weighted(ow),x = weighted(xw),y = weighted(yw),nx = x.norm
            let z = x.cross(y),nz = z.norm
            guard nx.isFinite,nz.isFinite,nx > 1e-12,nz > 1e-12*max(1,nx*y.norm) else {
                throw VivoChemistryError.invalid("degenerate local-coordinate virtual-site frame")
            }
            let ex = x/nx,ez = z/nz,ey = ez.cross(ex)
            position = origin+ex*local.x+ey*local.y+ez*local.z
            for i in p.indices { for axis in 0..<3 {
                let dx = axes[axis]*xw[i],dy = axes[axis]*yw[i]
                let dex = (dx-ex*ex.dot(dx))/nx
                let dz = dx.cross(y)+x.cross(dy),dez = (dz-ez*ez.dot(dz))/nz
                let dey = dez.cross(ex)+ez.cross(dex)
                columns[i][axis] = axes[axis]*ow[i]+dex*local.x+dey*local.y+dez*local.z
            } }
        }
        guard position.isFinite,columns.joined().allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("dependent-site construction overflow") }
        return (position,columns.map { .init(x: $0[0],y: $0[1],z: $0[2]) })
    }
}

public extension VivoClassicalSystem {
    /// Shared definition expansion. Full system validation is deliberately outside
    /// this helper so the canonical validator can validate the graph without recursion.
    func resolvedVirtualSiteGraph() throws -> VivoVirtualSiteGraph {
        guard !particles.contains(where: { $0.role == .drude }) else {
            throw VivoChemistryError.unsupported("Drude oscillators are not dependent virtual sites")
        }
        let graph = try VivoVirtualSiteGraph(particleCount: particles.count,
            physicalParticles: particles.filter { $0.role == .atom }.map(\.index),
            definitions: (linearVirtualSites ?? []).map(VivoDependentSite.init)+(virtualSiteDefinitions ?? []))
        guard Set(graph.sites.map(\.siteParticle)) == Set(particles.filter { $0.role == .virtualSite }.map(\.index)) else {
            throw VivoChemistryError.invalid("dependent-site role/definition ownership")
        }
        return graph
    }
}
