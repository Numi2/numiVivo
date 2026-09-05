import Foundation
#if NUMIVIVO_SCOPED_NUMERICS
@testable import NumiVivoNumerics
#endif

@main struct SmoothSolventChecks {
    static func main() throws {
        guard CommandLine.arguments.count==2 else{throw VivoChemistryError.invalid("output directory required")}
        let root=URL(fileURLWithPath:CommandLine.arguments[1]);try FileManager.default.createDirectory(at:root,withIntermediateDirectories:true)
        var count=0,values:[String:Double]=[:]
        func check(_ flag:Bool,_ label:String) throws{guard flag else{throw VivoChemistryError.invalid("solvent regression: \(label)")};count+=1;print("PASS \(label)")}
        for n in [50,110,302] {
            let grid=try VivoSolventLebedev.grid(n)
            let norm=grid.reduce(0.0){$0+$1.1},x2=grid.reduce(0.0){$0+$1.1*$1.0.x*$1.0.x}
            let x4=grid.reduce(0.0){$0+$1.1*pow($1.0.x,4)},xyz=grid.reduce(0.0){$0+$1.1*$1.0.x*$1.0.y*$1.0.z}
            try check(abs(norm-1)<1e-13 && abs(x2-1.0/3)<1e-13 && abs(x4-0.2)<1e-13 && abs(xyz)<1e-13,"Lebedev \(n) low-order spherical moments")
        }
        let system=VivoElectronicSystem(nuclei:[.init(atomicNumber:1,positionBohr:.init(0,0,-0.7)),.init(atomicNumber:1,positionBohr:.init(0,0,0.7))],alphaElectrons:1,betaElectrons:1)
        let basis=VivoGaussianBasis.hydrogenSTO3G(nucleusIndices:[0,1]),ao=try VivoGaussianIntegralEngine.compute(system:system,basis:basis)
        let hf=try VivoHartreeFock.solve(system:system,integrals:ao),op=try VivoSmoothCPCMOperator(system:system,basis:basis)
        var rotation=VivoQMMatrix(2,2);rotation[0,1]=0.24;rotation[1,0] = -0.24
        let c=try hf.alphaCoefficients.multiplied(by:VivoQMDenseAlgebra.orbitalRotation(generator:rotation))
        let density=VivoHartreeFock.density(c,occupied:1).scaled(2),result=try op.evaluate(totalDensity:density)
        var derivative=VivoQMMatrix(2,2)
        for p in 0..<2{for q in 0..<2{derivative[p,q]=2*(c[p,1]*c[q,0]+c[p,0]*c[q,1])}}
        let analytic=zip(result.reactionPotentialMatrix.values,derivative.values).reduce(0.0){$0+$1.0*$1.1}
        func energy(_ step:Double) throws -> Double {
            var k=VivoQMMatrix(2,2);k[0,1]=step;k[1,0] = -step
            let rotated=try c.multiplied(by:VivoQMDenseAlgebra.orbitalRotation(generator:k))
            return try op.evaluate(totalDensity:VivoHartreeFock.density(rotated,occupied:1).scaled(2)).polarizationEnergyHartree
        }
        let numeric=try (energy(1e-5)-energy(-1e-5))/2e-5
        values["variationalDerivativeError"]=abs(analytic-numeric)
        try check(abs(analytic-numeric)<1e-8,"reaction potential is the variational derivative of polarization energy")
        let offset=SIMD3<Double>(1.4,-2.1,0.33)
        let shifted=VivoElectronicSystem(nuclei:system.nuclei.map{.init(atomicNumber:$0.atomicNumber,positionBohr:$0.positionBohr+offset)},alphaElectrons:1,betaElectrons:1)
        let shiftedOp=try VivoSmoothCPCMOperator(system:shifted,basis:basis),translated=try shiftedOp.evaluate(totalDensity:density)
        values["translationEnergyError"]=abs(translated.polarizationEnergyHartree-result.polarizationEnergyHartree)
        try check(abs(translated.polarizationEnergyHartree-result.polarizationEnergyHartree)<1e-10,"translated cavity and solute retain energy")
        var badSystem=system;badSystem.pointCharges=[.init(chargeE:0.2,positionBohr:.init(100,0,0))]
        var rejected=false;do{_ = try VivoSmoothCPCMOperator(system:badSystem,basis:basis)}catch is VivoChemistryError{rejected=true}
        try check(rejected,"outlying point charge rejected without a declared enclosing cavity")
        var configuration=VivoSmoothCPCMConfiguration();configuration.maximumTesserae=1
        rejected=false;do{_ = try VivoSmoothCPCMOperator(system:system,basis:basis,configuration:configuration)}catch is VivoChemistryError{rejected=true}
        try check(rejected,"surface allocation is bounded")
        var wrong=density;wrong[0,0]+=0.1
        rejected=false;do{_ = try op.evaluate(totalDensity:wrong)}catch is VivoChemistryError{rejected=true}
        try check(rejected,"AO-metric electron-count mismatch rejected")
        values["trueRelativeResidual"]=result.trueRelativeResidual
        try check(result.trueRelativeResidual<5e-10,"original surface equation residual satisfies bound")
        struct Report:Codable{let checksPassed:Int;let observations:[String:Double];let scope:String}
        let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys,.prettyPrinted]
        try encoder.encode(Report(checksPassed:count,observations:values,scope:"smooth C-PCM electrostatic invariants and diagnostics; no thermal or nonelectrostatic-solvent model")).write(to:root.appendingPathComponent("results.json"),options:.atomic)
    }
}
