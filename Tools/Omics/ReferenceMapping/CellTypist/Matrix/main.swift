import Foundation
struct Row: Decodable { let indices: [Int]; let counts: [UInt64] }
struct Input: Decodable { let features: [String]; let rows: [Row] }
let model = try JSONDecoder().decode(VivoCellTypistReferenceModel.self, from: Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[1])))
let input = try JSONDecoder().decode(Input.self, from: Data(contentsOf: URL(fileURLWithPath: CommandLine.arguments[2])))
let predictor = try VivoCellTypistReference(model: model, sourceFeatures: input.features)
let totals = input.rows.map { $0.counts.reduce(UInt64(0), +) }
var entries: [(Int, Int, UInt64)] = []
for (row, values) in input.rows.enumerated() {
    for (feature, count) in zip(values.indices, values.counts) { entries.append((row, feature, count)) }
}
var maximumError = 0.0
for layout in ["csr", "csc"] {
    if layout == "csc" { entries.sort { ($0.1, $0.0) < ($1.1, $1.0) } }
    var seen = 0
    try predictor.predictSparseMatrix(rowTotals: totals, scan: { accept in
        for e in entries { try accept(e.0, e.1, e.2) }
    }, emit: { row, result in
        let expected = try predictor.predict(indices: input.rows[row].indices, counts: input.rows[row].counts)
        precondition(row == seen && result.label == expected.label)
        for (a, b) in zip(result.decisions, expected.decisions) { maximumError = max(maximumError, abs(a - b)) }
        seen += 1
    })
    precondition(seen == input.rows.count)
}
precondition(maximumError < 1e-10)
var rejected = 0
for kind in 0..<5 {
    var emitted = 0
    do {
        var badTotals = totals
        if kind == 0 { badTotals[0] += 1 }
        try predictor.predictSparseMatrix(rowTotals: badTotals, scan: { accept in
            if kind == 1 { try accept(-1, 0, 1) }
            if kind == 2 { try accept(0, input.features.count, 1) }
            if kind == 3 { try accept(0, 0, 0) }
            if kind == 4 { try accept(0, 0, 1); try accept(0, 0, 1) }
            for e in entries { try accept(e.0, e.1, e.2) }
        }, emit: { _, _ in emitted += 1 })
        fatalError("accepted invalid matrix")
    } catch { precondition(emitted == 0); rejected += 1 }
}
print("{\"cells\":\(totals.count),\"records\":\(entries.count),\"layouts\":[\"CSR\",\"CSC\"],\"maximumDecisionError\":\(maximumError),\"rejections\":\(rejected)}")
