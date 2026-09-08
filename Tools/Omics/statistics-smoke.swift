import Foundation
@main struct Smoke {
    static func main() throws {
        let x = [[1.0, 0], [1, 0], [1, 0], [1, 1], [1, 1], [1, 1]]
        let fit = try VivoOmicsQR(design: x).fit([1, 2, 3, 4, 5, 6], contrast: [0, 1])
        precondition(abs(fit.effect - 3) < 1e-12 && abs(fit.residualVariance - 1) < 1e-12)
        precondition(abs(fit.contrastVarianceScale - 2.0 / 3) < 1e-12)
        let p = try VivoOmicsLinearStatistics.studentTwoSidedP(t: 2.7764451051977987, degreesOfFreedom: 4)
        precondition(abs(p - 0.05) < 1e-10)
        let q = try VivoOmicsLinearStatistics.studentCriticalValue(degreesOfFreedom: 4)
        precondition(abs(q - 2.7764451051977987) < 1e-9)
        let adjusted = try VivoOmicsLinearStatistics.benjaminiHochberg([0.01, 0.04, 0.03, 0.8])
        precondition(zip(adjusted, [0.04, 0.05333333333333334, 0.05333333333333334, 0.8]).allSatisfy { abs($0 - $1) < 1e-14 })
        do { _ = try VivoOmicsQR(design: [[1, 0, 0], [1, 1, 1], [1, 0, 0], [1, 1, 1]]); fatalError("confounding accepted") }
        catch VivoOmicsStatisticsError.invalid(_) {}
        let prior = try VivoOmicsLinearStatistics.variancePrior((1...40).map { Double($0) / 20 }, residualDF: 8)
        precondition(prior.variance > 0 && prior.degreesOfFreedom > 0)
        for df in [1.0, 2, 3, 8, 30, 1000] {
            let v = try VivoOmicsLinearStatistics.studentCriticalValue(degreesOfFreedom: df)
            print("quantile", df, v)
        }
        print("Native omics statistics smoke passed")
    }
}
