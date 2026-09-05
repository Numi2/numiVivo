import Foundation

/// Native import of the molecular Gaussian subset of Basis Set Exchange JSON.
/// ECPs and spherical d/f shells fail explicitly; no basis representation is
/// silently substituted. Supplied data remains a separately fingerprinted input.
public enum VivoBasisSetExchange {
    private struct Document: Decodable {
        let name: String?
        let version: String?
        let elements: [String:Element]
    }
    private struct Element: Decodable {
        let electron_shells: [Shell]?
        let ecp_electrons: Int?
    }
    private struct Shell: Decodable {
        let function_type: String
        let angular_momentum: [Int]
        let exponents: [String]
        let coefficients: [[String]]
    }
    public static func decode(_ data: Data, for system: VivoElectronicSystem, identifier: String,
                              source: String, budget: VivoChemistryBudget = .init()) throws -> VivoGaussianBasis {
        try budget.validate();try system.validate()
        guard data.count<=budget.maximumBytes,!identifier.isEmpty,!source.isEmpty else { throw VivoChemistryError.invalid("basis import size or provenance") }
        let document=try JSONDecoder().decode(Document.self,from:data)
        var shells:[VivoGaussianShell]=[]
        for (index,nucleus) in system.nuclei.enumerated() {
            guard let element=document.elements[String(nucleus.atomicNumber)],let definitions=element.electron_shells else {
                throw VivoChemistryError.invalid("basis has no all-electron shells for Z=\(nucleus.atomicNumber)")
            }
            guard (element.ecp_electrons ?? 0)==0 else { throw VivoChemistryError.unsupported("effective core potentials") }
            for definition in definitions {
                let angular=definition.angular_momentum
                guard !angular.isEmpty,angular.allSatisfy({(0...3).contains($0)}),Set(angular).count==angular.count,
                      !definition.exponents.isEmpty,definition.exponents.count<=128,!definition.coefficients.isEmpty else {
                    throw VivoChemistryError.invalid("BSE shell shape/angular momenta")
                }
                guard ["gto","gto_cartesian","gto_spherical"].contains(definition.function_type) else {
                    throw VivoChemistryError.unsupported("non-Gaussian BSE shell type")
                }
                if definition.function_type=="gto_spherical" && angular.contains(where:{$0>=2}) {
                    throw VivoChemistryError.unsupported("spherical d/f shell transformation is not implemented")
                }
                guard angular.count==1 || angular.count==definition.coefficients.count else {
                    throw VivoChemistryError.unsupported("ambiguous combined-shell general contractions")
                }
                let exponents=try definition.exponents.map { text -> Double in
                    guard let value=Double(text.replacingOccurrences(of:"D",with:"E").replacingOccurrences(of:"d",with:"e")),value.isFinite,value>0 else {
                        throw VivoChemistryError.invalid("BSE Gaussian exponent")
                    };return value
                }
                for (column,coefficients) in definition.coefficients.enumerated() {
                    guard coefficients.count==exponents.count else { throw VivoChemistryError.invalid("BSE contraction length") }
                    let values=try coefficients.map { text -> Double in
                        guard let value=Double(text.replacingOccurrences(of:"D",with:"E").replacingOccurrences(of:"d",with:"e")),value.isFinite else {
                            throw VivoChemistryError.invalid("BSE contraction coefficient")
                        };return value
                    }
                    shells.append(.init(nucleusIndex:index,angularMomentum:angular.count==1 ? angular[0] : angular[column],
                                        primitives:zip(exponents,values).map{.init(exponent:$0.0,coefficient:$0.1)}))
                }
            }
        }
        let basis=VivoGaussianBasis(identifier:identifier,shells:shells,
            source:source+"; BSE name="+(document.name ?? "unspecified")+"; version="+(document.version ?? "unspecified"))
        _ = try VivoGaussianIntegralEngine.expanded(system:system,basis:basis,budget:budget)
        return basis
    }
}
