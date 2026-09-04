import Foundation
@main struct CouplingIterationChecks {
    static func main() throws {
        var checked = 0
        func expect(_ value: @autoclosure () throws -> Bool) throws {
            guard try value() else { throw NSError(domain: "coupling-check", code: 1) }
            checked += 1
        }
        func rejects(_ action: () throws -> Void) throws {
            do { try action() } catch { checked += 1; return }
            throw NSError(domain: "coupling-check-expected-rejection", code: 2)
        }
        try expect(VivoCouplingIterationController.residual(candidate: [], reference: []) == 0)
        try expect(VivoCouplingIterationController.residual(candidate: [0,0], reference: [0,0]) == 0)
        try expect(VivoCouplingIterationController.residual(candidate: [0], reference: [1]) == 1)
        let a = try VivoCouplingIterationController.residual(candidate: [2,8], reference: [1,4])
        let b = try VivoCouplingIterationController.residual(candidate: [2000,0.008], reference: [1000,0.004])
        try expect(abs(a-b) < 1e-6)
        try rejects { _ = try VivoCouplingIterationController.residual(candidate: [1], reference: []) }
        try rejects { _ = try VivoCouplingIterationController.residual(candidate: [.nan], reference: [1]) }
        try rejects { _ = try VivoCouplingIterationController(relaxation: 0) }
        var controller = try VivoCouplingIterationController(relaxation: 0.7)
        var value: [Float] = [0]
        value = try controller.next(current: value, proposed: [10])
        try expect(abs(value[0] - 7) < 1e-6)
        value = try controller.next(current: value, proposed: [10])
        try expect(abs(value[0] - 10) < 1e-5)
        try expect(VivoCouplingIterationController.residual(candidate: [10], reference: value) < 1e-5)
        var fixed = try VivoCouplingIterationController(relaxation: 0.7)
        var x: [Float] = [0]
        for _ in 0..<30 { x = try fixed.next(current: x, proposed: [0.5 * x[0] + 2]) }
        try expect(abs(x[0] - 4) < 1e-5)
        try rejects { _ = try fixed.next(current: [0], proposed: [.infinity]) }
        print("Coupling iteration portable checks passed: \(checked)")
    }
}
