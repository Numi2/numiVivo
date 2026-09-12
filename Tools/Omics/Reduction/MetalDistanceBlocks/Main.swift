import Foundation
import Metal
let begin=Date()
let args=CommandLine.arguments
precondition(args.count==6)
let mode=args[1], source=URL(fileURLWithPath:args[2]), n=Int(args[3])!, d=Int(args[4])!, output=args[5]
precondition((20...50000).contains(n) && (1...64).contains(d) && ["cpu","metal"].contains(mode))
precondition(!FileManager.default.fileExists(atPath:output))
let data=try Data(contentsOf:source,options:.mappedIfSafe)
precondition(data.count==n*d*16)
var scores=[Float](repeating:0,count:n*d)
data.withUnsafeBytes { b in
 for i in 0..<n*d {
  precondition(b.loadUnaligned(fromByteOffset:i*16,as:UInt32.self)==UInt32(i/d))
  precondition(b.loadUnaligned(fromByteOffset:i*16+4,as:UInt32.self)==UInt32(i%d))
  let x=b.loadUnaligned(fromByteOffset:i*16+8,as:Double.self);precondition(x.isFinite)
  scores[i]=Float(x);precondition(scores[i].isFinite)
 }
}
let shader="""
#include <metal_stdlib>
using namespace metal;
kernel void distances(device const float *x [[buffer(0)]], device float *y [[buffer(1)]], constant uint4 &a [[buffer(2)]], uint2 id [[thread_position_in_grid]]) {
 if(id.x>=a.x || id.y>=a.w) return;
 float s=0;
 for(uint c=0;c<a.y;c++) { float z=x[(a.z+id.y)*a.y+c]-x[id.x*a.y+c];s=s+z*z; }
 y[id.y*a.x+id.x]=s;
}
"""
let tile=64, k=20
var device:MTLDevice?, queue:MTLCommandQueue?, pipeline:MTLComputePipelineState?, input:MTLBuffer?, distances:MTLBuffer?
if mode=="metal" {
 device=MTLCreateSystemDefaultDevice();precondition(device != nil && device!.hasUnifiedMemory && !device!.name.lowercased().contains("paravirtual"))
 queue=device!.makeCommandQueue();let options=MTLCompileOptions();options.fastMathEnabled=false
 let library=try device!.makeLibrary(source:shader,options:options)
 pipeline=try device!.makeComputePipelineState(function:library.makeFunction(name:"distances")!)
 input=scores.withUnsafeBytes {device!.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)}
 distances=device!.makeBuffer(length:tile*n*4,options:.storageModeShared)
}
var ids=[UInt32]();var values=[Float]();ids.reserveCapacity(n*k);values.reserveCapacity(n*k)
var cpu=[Float](repeating:0,count:mode=="cpu" ? tile*n : 0)
var gpuSeconds=0.0
for start in stride(from:0,to:n,by:tile) {
 let count=min(tile,n-start)
 if mode=="metal" {
  let command=queue!.makeCommandBuffer()!, encoder=command.makeComputeCommandEncoder()!
  encoder.setComputePipelineState(pipeline!);encoder.setBuffer(input,offset:0,index:0);encoder.setBuffer(distances,offset:0,index:1)
  var axes=SIMD4<UInt32>(UInt32(n),UInt32(d),UInt32(start),UInt32(count));encoder.setBytes(&axes,length:16,index:2)
  encoder.dispatchThreads(MTLSize(width:n,height:count,depth:1),threadsPerThreadgroup:MTLSize(width:32,height:4,depth:1))
  encoder.endEncoding();command.commit();command.waitUntilCompleted();precondition(command.status == .completed)
  gpuSeconds += command.gpuEndTime-command.gpuStartTime
 } else {
  for row in 0..<count {for other in 0..<n {
   var sum:Float=0
   for c in 0..<d {let delta=scores[(start+row)*d+c]-scores[other*d+c];sum += delta*delta}
   cpu[row*n+other]=sum
  }}
 }
 let metalOutput = distances?.contents().assumingMemoryBound(to:Float.self)
 for row in 0..<count {
  var best=[(Float,UInt32)]();best.reserveCapacity(k-1)
  for other in 0..<n where other != start+row {
   let x:Float=mode=="metal" ? metalOutput![row*n+other] : cpu[row*n+other]
   precondition(x.isFinite && x>=0)
   let id=UInt32(other)
   if best.count==k-1 {let last=best[k-2];if x>last.0 || (x==last.0 && id>last.1) {continue}}
   let pos=best.firstIndex(where:{ x<$0.0 || (x==$0.0 && id<$0.1) }) ?? best.count
   best.insert((x,id),at:pos);if best.count>k-1 {best.removeLast()}
  }
  ids.append(UInt32(start+row));values.append(0)
  for (x,id) in best {ids.append(id);values.append(x)}
 }
}
try FileManager.default.createDirectory(atPath:output,withIntermediateDirectories:false)
try ids.withUnsafeBytes { try Data($0).write(to:URL(fileURLWithPath:output+"/indices.bin")) }
try values.withUnsafeBytes {try Data($0).write(to:URL(fileURLWithPath:output+"/squared-distances.bin"))}
let report:[String:Any]=["backend":mode,"cells":n,"dimensions":d,"neighborsIncludingSelf":k,"queryTile":tile,"distanceBufferBytes":tile*n*4,"secondsIncludingLoadCompileComputeSelectionWrite":Date().timeIntervalSince(begin),"gpuCommandSeconds":gpuSeconds,"device":device?.name ?? "CPU","scope":"Research FP32 kNN distance blocks; not native product integration or biological qualification"]
let reportData=try JSONSerialization.data(withJSONObject:report,options:[.sortedKeys,.prettyPrinted]);try reportData.write(to:URL(fileURLWithPath:output+"/report.json"));print(String(data:reportData,encoding:.utf8)!)
