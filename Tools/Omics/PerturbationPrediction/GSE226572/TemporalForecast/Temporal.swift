import Foundation
struct Observation:Codable {let donor:String,time:Double,response:[Double]}
struct Input:Codable {let donor:String,time:Double,control:[Double],training:[Observation]}
struct Output:Codable {let predicted:[Double],baseline:[Double],lowerTime:Double,upperTime:Double}
@main struct Main {
 static func main() throws {
  let input=try JSONDecoder().decode(Input.self,from:Data(contentsOf:URL(fileURLWithPath:CommandLine.arguments[1])))
  let g=input.control.count
  guard g>0,g<=100000,input.time>0,input.time.isFinite,input.training.count<=256,
   input.control.allSatisfy({$0.isFinite && $0>=0}),
   input.training.allSatisfy({$0.donor != input.donor && $0.time < input.time && $0.time>0 && $0.time.isFinite && $0.response.count==g && $0.response.allSatisfy(\.isFinite)}) else {throw NSError(domain:"invalid or leaking temporal input",code:1)}
  if input.training.isEmpty {
   let e=JSONEncoder();e.outputFormatting=[.sortedKeys]
   try e.encode(Output(predicted:input.control,baseline:input.control,lowerTime:0,upperTime:0)).write(to:URL(fileURLWithPath:CommandLine.arguments[2]),options:.withoutOverwriting)
   return
  }
  let times=Set(input.training.map(\.time)).sorted(),donors=Set(input.training.map(\.donor)).sorted()
  var curve=[Double:[Double]]();curve[0]=[Double](repeating:0,count:g)
  for time in times {
   let obs=input.training.filter{$0.time==time};guard Set(obs.map(\.donor)).count==obs.count else {throw NSError(domain:"duplicate donor time",code:2)}
   var mean=[Double](repeating:0,count:g)
   for row in obs {for j in 0..<g {mean[j]+=row.response[j]/Double(obs.count)}}
   curve[time]=mean
  }
  let knots=[0.0]+times,lo=knots.last(where:{$0<=input.time}) ?? 0,hi=knots.first(where:{$0>=input.time}) ?? times.last!
  let lower=min(lo,hi),upper=hi,alpha=upper==lower ? 0 : (input.time-lower)/(upper-lower)
  var predicted=input.control,baseline=input.control,mean=[Double](repeating:0,count:g)
  for donor in donors {
   let obs=input.training.filter{$0.donor==donor}
   for row in obs {for j in 0..<g {mean[j]+=row.response[j]/Double(donors.count*obs.count)}}
  }
  for j in 0..<g {
   predicted[j]=max(0,input.control[j]+(1-alpha)*curve[lower]![j]+alpha*curve[upper]![j])
   baseline[j]=max(0,input.control[j]+mean[j])
  }
  let e=JSONEncoder();e.outputFormatting=[.sortedKeys]
  try e.encode(Output(predicted:predicted,baseline:baseline,lowerTime:lower,upperTime:upper)).write(to:URL(fileURLWithPath:CommandLine.arguments[2]),options:.withoutOverwriting)
 }
}
