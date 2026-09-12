import Foundation
@main struct Tests {
 static func main() throws {
  let ids=["g"], donors=["a","b","c"]
  let m=try VivoResponseShrinkage.fit(featureIDs:ids,trainingDonorIDs:donors,controls:[[0],[1],[2]],treated:[[2],[3],[4]])
  precondition(m.selectedPenalty == nil && m.responseSlope == [0])
  precondition(tryValue(m,3)==5)
  let linear=try VivoResponseShrinkage.fit(featureIDs:ids,trainingDonorIDs:donors,controls:[[0],[2],[4]],treated:[[1],[4],[7]],penalties:[0])
  precondition(abs(tryValue(linear,6)-10)<1e-12)
  let clip=try VivoResponseShrinkage.fit(featureIDs:ids,trainingDonorIDs:donors,controls:[[2],[3],[4]],treated:[[0],[1],[2]])
  precondition(tryValue(clip,0)==0)
  let constant=try VivoResponseShrinkage.fit(featureIDs:ids,trainingDonorIDs:donors,controls:[[1],[1],[1]],treated:[[1],[2],[3]],penalties:[0])
  precondition(constant.responseSlope == [0])
  var rejected=0
  let cases:[() throws -> Void]=[
   { _=try VivoResponseShrinkage.predict(m,featureIDs:["other"],queryDonorID:"q",control:[1]) },
   { _=try VivoResponseShrinkage.predict(m,featureIDs:ids,queryDonorID:"a",control:[1]) },
   { _=try VivoResponseShrinkage.predict(m,featureIDs:ids,queryDonorID:"q",control:[.nan]) },
   { _=try VivoResponseShrinkage.fit(featureIDs:ids,trainingDonorIDs:["a","a","c"],controls:[[0],[1],[2]],treated:[[2],[3],[4]]) },
   { _=try VivoResponseShrinkage.fit(featureIDs:ids,trainingDonorIDs:donors,controls:[[0],[1],[2]],treated:[[2],[3],[4]],penalties:[1,0]) },
   { _=try VivoResponseShrinkage.fit(featureIDs:ids,trainingDonorIDs:donors,controls:[[0],[1],[2]],treated:[[2],[-1],[4]]) },
   { _=try VivoResponseShrinkage.fit(featureIDs:ids,trainingDonorIDs:["a","b"],controls:[[0],[1]],treated:[[2],[3]]) }
  ]
  for action in cases { do { try action() } catch { rejected+=1 } }
  precondition(rejected==cases.count);print("PASS: four analytical cases and seven input rejections")
 }
 static func tryValue(_ model: VivoResponseShrinkage.Model,_ value:Double)->Double { try! VivoResponseShrinkage.predict(model,featureIDs:["g"],queryDonorID:"q",control:[value])[0] }
}
