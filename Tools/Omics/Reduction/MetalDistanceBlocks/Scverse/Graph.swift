import Foundation
@testable import NumiVivoKit
let a=CommandLine.arguments;precondition(a.count==5)
let mode=a[1],source=URL(fileURLWithPath:a[2]),out=URL(fileURLWithPath:a[4])
let meta=try JSONDecoder().decode(VivoSingleCellCountMetadata.self,from:Data(contentsOf:URL(fileURLWithPath:a[3])))
let cells=meta.cells.map{VivoOmicsCellIdentity(sampleID:$0.sampleID,barcode:$0.barcode)}
var options=VivoSingleCellNeighborOptions();options.neighbors=20;options.maximumDistancePairs=500_000_000;options.representation = .pca
var execution=VivoPCANeighborExecution();execution.workers=1;execution.queryBlockRows=128;execution.candidateBlockRows=8192
precondition(mode=="cpu" || mode=="metal");if mode=="metal" {execution.backend = .metalFP32}
let start=Date()
let graph=try VivoWindowedPCANeighbors.run(source:source,cells:cells,dimensions:20,options:options,execution:execution)
let elapsed=Date().timeIntervalSince(start)
try FileManager.default.createDirectory(at:out,withIntermediateDirectories:false)
let bytes=try VivoCanonicalJSON.encode(graph);try bytes.write(to:out.appendingPathComponent("graph.json"))
let report=try JSONSerialization.data(withJSONObject:["backend":mode,"graphStageSeconds":elapsed,"cells":cells.count,"sourceScoreReadsIncluded":true],options:.sortedKeys)
try report.write(to:out.appendingPathComponent("timing.json"));print(String(data:report,encoding:.utf8)!)
