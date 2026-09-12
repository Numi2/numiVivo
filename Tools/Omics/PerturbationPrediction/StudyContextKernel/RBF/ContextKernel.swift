import Foundation
struct Input: Codable {
 let featureIDs:[String], trainingDonorIDs:[String], trainingStudies:[String]
 let controls:[[Double]], treated:[[Double]]
 let queryStudy:String, queryDonorIDs:[String], queryControls:[[Double]]
}
struct Model: Codable {
 let penalties:[Double?],selectedPenalty:Double?,innerStudyMSE:[[Double]]
 let trainingStudies:[String],trainingDonorIDs:[String],donorWeights:[Double]
 let queryResponseWeights:[[Double]],maximumSolveResidual:Double
}
struct Output: Codable {
 let featureIDs:[String],queryStudy:String,queryDonorIDs:[String],model:Model
 let predictedTreated:[[Double]],trainingMeanPredictions:[[Double]]
}
enum Failure:Error {case invalid(String)}
func weights(_ studies:[String])->[Double] {
 let unique=Set(studies);return studies.map {s in 1/Double(unique.count*studies.filter{$0==s}.count)}
}
struct Kernel {
 let w:[Double],center:[Double],x:[[Double]],gram:[[Double]]
 let bandwidthSquared:Double,rowMeans:[Double],grandMean:Double
 init(_ controls:[[Double]],_ studies:[String]) {
  let n=controls.count,g=controls[0].count;let weights=weights(studies);w=weights
  center=[Double](repeating:0,count:g);x=controls
  var d=Array(repeating:Array(repeating:0.0,count:n),count:n),pairs=[Double]()
  for i in 0..<n {for j in 0..<i {
   var distance=0.0;for c in 0..<g {let v=controls[i][c]-controls[j][c];distance+=v*v/Double(g)}
   d[i][j]=distance;d[j][i]=distance;pairs.append(distance)
  }}
  pairs.sort();let median=pairs.count%2==1 ? pairs[pairs.count/2] : (pairs[pairs.count/2-1]+pairs[pairs.count/2])/2
  let bandwidth=median>0 ? median : 1.0;bandwidthSquared=bandwidth
  let raw=d.map { $0.map {exp(-$0/(2*bandwidth))} }
  let means=raw.map {row in zip(row,weights).reduce(0.0){$0+$1.0*$1.1}}
  let grand=zip(means,weights).reduce(0.0){$0+$1.0*$1.1};rowMeans=means;grandMean=grand
  var k=raw
  for i in 0..<n {for j in 0..<n {k[i][j]=sqrt(weights[i]*weights[j])*(raw[i][j]-means[i]-means[j]+grand)}}
  gram=k
 }
 func responseWeights(_ queries:[[Double]],_ penalty:Double?) throws -> ([[Double]],Double) {
  guard let penalty else {return (queries.map {_ in w},0)}
  let n=w.count,g=center.count;var lower=gram
  for i in 0..<n {lower[i][i]+=penalty}
  for i in 0..<n {for j in 0...i {
   var value=lower[i][j];for k in 0..<j {value-=lower[i][k]*lower[j][k]}
   if i==j {guard value>0 && value.isFinite else {throw Failure.invalid("kernel Cholesky")};lower[i][j]=sqrt(value)}
   else {lower[i][j]=value/lower[j][j]}
  }}
  var out=[[Double]](),maxResidual=0.0
  for q in queries {
   try Task.checkCancellation();var rhs=[Double](repeating:0,count:n)
   var raw=[Double](repeating:0,count:n)
   for i in 0..<n {var distance=0.0;for j in 0..<g {let v=x[i][j]-q[j];distance+=v*v/Double(g)};raw[i]=exp(-distance/(2*bandwidthSquared))}
   let mean=zip(raw,w).reduce(0.0){$0+$1.0*$1.1}
   for i in 0..<n {rhs[i]=sqrt(w[i])*(raw[i]-rowMeans[i]-mean+grandMean)}
   var a=rhs
   for i in 0..<n {for j in 0..<i {a[i]-=lower[i][j]*a[j]};a[i]/=lower[i][i]}
   for i in stride(from:n-1,through:0,by:-1) {for j in (i+1)..<n {a[i]-=lower[j][i]*a[j]};a[i]/=lower[i][i]}
   for i in 0..<n {var value=penalty*a[i];for j in 0..<n {value+=gram[i][j]*a[j]};maxResidual=max(maxResidual,abs(value-rhs[i])/(1+abs(rhs[i])))}
   let effect=zip(w,a).map {sqrt($0)*$1},sum=effect.reduce(0,+)
   let beta=w.indices.map {w[$0]+effect[$0]-sum*w[$0]}
   guard beta.allSatisfy(\.isFinite),abs(beta.reduce(0,+)-1)<1e-9 else {throw Failure.invalid("response weights")}
   out.append(beta)
  }
  guard maxResidual<1e-9 else {throw Failure.invalid("kernel residual")}
  return(out,maxResidual)
 }
}
func predictions(_ controls:[[Double]],_ treated:[[Double]],_ query:[[Double]],_ coefficients:[[Double]]) throws -> [[Double]] {
 var out=[[Double]]();let n=controls.count,g=controls[0].count
 for i in query.indices {var row=query[i];for j in 0..<g {
  var response=0.0;for k in 0..<n {response+=coefficients[i][k]*(treated[k][j]-controls[k][j])}
  let value=row[j]+response;guard value.isFinite else {throw Failure.invalid("prediction")};row[j]=max(0,value)
 };out.append(row)}
 return out
}
func run(_ a:Input)throws->Output {
 let n=a.trainingDonorIDs.count,g=a.featureIDs.count,studies=Set(a.trainingStudies).sorted(),p:[Double?]=[0.01,0.1,1,10,100,nil]
 guard n>=4,n<=128,g>0,g<=100_000,n*g<=2_000_000,studies.count>=2,a.trainingStudies.count==n,!studies.contains(a.queryStudy),
 Set(a.trainingDonorIDs).count==n,Set(a.queryDonorIDs).count==a.queryDonorIDs.count,!a.queryDonorIDs.isEmpty,a.queryDonorIDs.count<=128,
 Set(a.trainingDonorIDs).isDisjoint(with:Set(a.queryDonorIDs)),Set(a.featureIDs).count==g,
 a.controls.count==n,a.treated.count==n,a.queryControls.count==a.queryDonorIDs.count,a.queryDonorIDs.count*g<=2_000_000,
 (a.controls+a.treated+a.queryControls).allSatisfy({$0.count==g && $0.allSatisfy{$0.isFinite && $0>=0 && $0<=100}}),
 (a.trainingDonorIDs+a.queryDonorIDs+a.trainingStudies+a.featureIDs+[a.queryStudy]).allSatisfy({!$0.isEmpty && $0.utf8.count<=1024})
 else {throw Failure.invalid("input identities, study exclusion, bounds or log CPM")}
 var losses=[[Double]](),maxResidual=0.0
 for study in studies {
  let keep=(0..<n).filter{a.trainingStudies[$0] != study},test=(0..<n).filter{a.trainingStudies[$0]==study}
  guard keep.count>=2 else {throw Failure.invalid("inner training donors")}
  let x=keep.map{a.controls[$0]},y=keep.map{a.treated[$0]},q=test.map{a.controls[$0]},kernel=Kernel(x,keep.map{a.trainingStudies[$0]})
  var loss=[Double]()
  for penalty in p {
   let (beta,residual)=try kernel.responseWeights(q,penalty);maxResidual=max(maxResidual,residual)
   let predicted=try predictions(x,y,q,beta);var error=0.0
   for i in test.indices {for j in 0..<g {let e=predicted[i][j]-a.treated[test[i]][j];error+=e*e/Double(test.count*g)}}
   guard error.isFinite else {throw Failure.invalid("loss")};loss.append(error)
  }
  losses.append(loss)
 }
 let means=p.indices.map{k in losses.reduce(0){$0+$1[k]}/Double(studies.count)}
 var best=0;for k in p.indices where means[k]<=means[best] {best=k}
 let kernel=Kernel(a.controls,a.trainingStudies), (beta,residual)=try kernel.responseWeights(a.queryControls,p[best]);maxResidual=max(maxResidual,residual)
 let predicted=try predictions(a.controls,a.treated,a.queryControls,beta)
 let baseline=try predictions(a.controls,a.treated,a.queryControls,a.queryControls.map{_ in kernel.w})
 let model=Model(penalties:p,selectedPenalty:p[best],innerStudyMSE:losses,trainingStudies:studies,trainingDonorIDs:a.trainingDonorIDs,donorWeights:kernel.w,queryResponseWeights:beta,maximumSolveResidual:maxResidual)
 return Output(featureIDs:a.featureIDs,queryStudy:a.queryStudy,queryDonorIDs:a.queryDonorIDs,model:model,predictedTreated:predicted,trainingMeanPredictions:baseline)
}
@main struct Driver {static func main(){do {
 let bytes=try Data(contentsOf:URL(fileURLWithPath:CommandLine.arguments[1]))
 guard let object=try JSONSerialization.jsonObject(with:bytes) as? [String:Any],Set(object.keys)==Set(["featureIDs","trainingDonorIDs","trainingStudies","controls","treated","queryStudy","queryDonorIDs","queryControls"]) else {throw Failure.invalid("input fields; query outcomes forbidden")}
 let a=try JSONDecoder().decode(Input.self,from:bytes),b=try run(a);let e=JSONEncoder();e.outputFormatting=[.sortedKeys];try e.encode(b).write(to:URL(fileURLWithPath:CommandLine.arguments[2]),options:.withoutOverwriting)
} catch {FileHandle.standardError.write(Data("\(error)\n".utf8));exit(2)}}}
