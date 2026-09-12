import Foundation
struct Input: Codable {
 let featureIDs:[String], trainingDonorIDs:[String], trainingStudies:[String]
 let controls:[[Double]], treated:[[Double]]
 let queryStudy:String, queryDonorIDs:[String], queryControls:[[Double]]
}
struct Model: Codable {
 let penalties:[Double?], selectedPenalty:Double?, innerStudyMSE:[[Double]]
 let trainingStudies:[String], trainingDonorIDs:[String], donorWeights:[Double]
 let controlMean:[Double], responseMean:[Double], responseSlope:[Double]
}
struct Output: Codable { let featureIDs:[String],queryStudy:String,queryDonorIDs:[String];let model:Model;let predictedTreated:[[Double]],trainingMeanPredictions:[[Double]] }
enum Failure:Error { case invalid(String) }
func weights(_ studies:[String])->[Double] {
 let unique=Set(studies);return studies.map { s in 1/Double(unique.count*studies.filter { $0==s }.count) }
}
func coefficients(_ x:[Double],_ d:[Double],_ w:[Double],_ lam:Double?)->(Double,Double,Double) {
 let center=zip(x,w).reduce(0) { $0+$1.0*$1.1 },mean=zip(d,w).reduce(0) { $0+$1.0*$1.1 }
 guard let lam else {return(center,mean,0)}
 var cov=0.0,variance=0.0
 for i in x.indices { let v=x[i]-center;cov+=w[i]*v*(d[i]-mean);variance+=w[i]*v*v }
 let den=variance+lam*(1-w.reduce(0) { $0+$1*$1 })
 return(center,mean,den>0 ? cov/den : 0)
}
func predict(_ q:Double,_ c:Double,_ d:Double,_ b:Double)throws->Double {
 let raw=q+d+b*(q-c);guard raw.isFinite else {throw Failure.invalid("nonfinite prediction")};return max(0,raw)
}
func run(_ a:Input)throws->Output {
 let n=a.trainingDonorIDs.count,g=a.featureIDs.count,studies=Set(a.trainingStudies).sorted(),p:[Double?]=[0,0.01,0.1,1,10,100,nil]
 guard n>=4,n<=128,g>0,g<=100_000,n*g<=2_000_000,studies.count>=2,a.trainingStudies.count==n,!studies.contains(a.queryStudy),
       Set(a.trainingDonorIDs).count==n,Set(a.queryDonorIDs).count==a.queryDonorIDs.count,!a.queryDonorIDs.isEmpty,
       Set(a.trainingDonorIDs).isDisjoint(with:Set(a.queryDonorIDs)),Set(a.featureIDs).count==g,
       a.controls.count==n,a.treated.count==n,a.queryControls.count==a.queryDonorIDs.count,a.queryDonorIDs.count*g<=2_000_000,
       (a.controls+a.treated+a.queryControls).allSatisfy({$0.count==g && $0.allSatisfy {$0.isFinite && $0>=0 && $0<=100}}),
       (a.trainingDonorIDs+a.queryDonorIDs+a.trainingStudies+a.featureIDs+[a.queryStudy]).allSatisfy({!$0.isEmpty && $0.utf8.count<=1024}) else {throw Failure.invalid("input identities, study exclusion, bounds or log CPM")}
 let splits=studies.map { s in ((0..<n).filter {a.trainingStudies[$0] != s},(0..<n).filter {a.trainingStudies[$0]==s}) }
 guard splits.allSatisfy({$0.0.count>=2 && !$0.1.isEmpty}) else {throw Failure.invalid("inner training donors")}
 var loss=Array(repeating:Array(repeating:0.0,count:p.count),count:studies.count)
 for gene in 0..<g {
  if gene%256==0 {try Task.checkCancellation()}
  let x=a.controls.map {$0[gene]},d=zip(a.controls,a.treated).map {$1[gene]-$0[gene]}
  for (si,split) in splits.enumerated() {
   let keep=split.0,test=split.1,xx=keep.map {x[$0]},dd=keep.map {d[$0]},ww=weights(keep.map {a.trainingStudies[$0]})
   for k in p.indices {
    let(c,m,b)=coefficients(xx,dd,ww,p[k])
    for j in test { let err=try predict(x[j],c,m,b)-a.treated[j][gene];loss[si][k]+=err*err/Double(test.count*g) }
   }
  }
 }
 let means=p.indices.map { k in loss.reduce(0) {$0+$1[k]}/Double(studies.count) }
 guard means.allSatisfy(\.isFinite) else {throw Failure.invalid("loss")}
 var best=0;for k in p.indices where means[k]<=means[best] {best=k}
 let w=weights(a.trainingStudies);var centers=[Double](),responses=[Double](),slopes=[Double]()
 for gene in 0..<g {let x=a.controls.map {$0[gene]},d=zip(a.controls,a.treated).map {$1[gene]-$0[gene]};let(c,m,b)=coefficients(x,d,w,p[best]);centers.append(c);responses.append(m);slopes.append(b)}
 var output=[[Double]](),baseline=[[Double]]()
 for q in a.queryControls {var o=[Double](),base=[Double]();for gene in 0..<g {o.append(try predict(q[gene],centers[gene],responses[gene],slopes[gene]));base.append(try predict(q[gene],centers[gene],responses[gene],0))};output.append(o);baseline.append(base)}
 let model=Model(penalties:p,selectedPenalty:p[best],innerStudyMSE:loss,trainingStudies:studies,trainingDonorIDs:a.trainingDonorIDs,donorWeights:w,controlMean:centers,responseMean:responses,responseSlope:slopes)
 return Output(featureIDs:a.featureIDs,queryStudy:a.queryStudy,queryDonorIDs:a.queryDonorIDs,model:model,predictedTreated:output,trainingMeanPredictions:baseline)
}
@main struct Driver {static func main()throws {let a=try JSONDecoder().decode(Input.self,from:Data(contentsOf:URL(fileURLWithPath:CommandLine.arguments[1])));let b=try run(a);let e=JSONEncoder();e.outputFormatting=[.sortedKeys];try e.encode(b).write(to:URL(fileURLWithPath:CommandLine.arguments[2]),options:.withoutOverwriting)}}
