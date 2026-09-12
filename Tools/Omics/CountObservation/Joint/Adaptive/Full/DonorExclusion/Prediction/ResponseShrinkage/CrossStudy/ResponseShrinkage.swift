import Foundation

/// Donor-pseudobulk response regression. Inputs are log1p CPM with full-library
/// denominators; no single-cell likelihood or uncertainty calibration is implied.
public enum VivoResponseShrinkage {
    public enum Failure: Error { case invalidInput(String) }
    public struct Model: Codable, Sendable {
        public let featureIDs: [String]
        public let trainingDonorIDs: [String]
        public let selectedPenalty: Double?
        public let penalties: [Double?]
        public let innerDonorMSE: [Double]
        public let controlMean: [Double]
        public let responseMean: [Double]
        public let responseSlope: [Double]
    }
    private static func validIDs(_ ids: [String]) -> Bool {
        !ids.isEmpty && Set(ids).count == ids.count && ids.allSatisfy { !$0.isEmpty && $0.utf8.count <= 1024 }
    }
    private static func sum(_ xs: [Double]) -> Double {
        var total=0.0, correction=0.0
        for x in xs { let y=x-correction; let next=total+y; correction=(next-total)-y; total=next }
        return total
    }
    private static func coefficients(_ x: [Double], _ d: [Double], _ penalty: Double?) -> (Double,Double,Double) {
        let xm=sum(x)/Double(x.count), dm=sum(d)/Double(d.count)
        guard let penalty else { return (xm,dm,0) }
        let xc=x.map { $0-xm }
        let denominator=sum(xc.map { $0*$0 })+penalty*Double(x.count-1)
        let numerator=sum(zip(xc,d).map { $0*($1-dm) })
        return (xm,dm,denominator>0 ? numerator/denominator : 0)
    }
    public static func fit(featureIDs: [String], trainingDonorIDs: [String], controls: [[Double]], treated: [[Double]], penalties: [Double?] = [0,0.01,0.1,1,10,100,nil]) throws -> Model {
        let n=trainingDonorIDs.count, g=featureIDs.count
        guard n>=3,n<=128,g<=100_000,n*g<=2_000_000,validIDs(featureIDs),validIDs(trainingDonorIDs),controls.count==n,treated.count==n,
              (controls+treated).allSatisfy({ $0.count==g && $0.allSatisfy { $0.isFinite && $0>=0 && $0<=100 } }),
              !penalties.isEmpty,penalties.count<=32,penalties.allSatisfy({ $0 == nil || ($0!.isFinite && $0!>=0 && $0!<=1e6) }) else { throw Failure.invalidInput("training dimensions, identities, log CPM or penalties") }
        // The candidate order is ascending regularization. Nil denotes infinity.
        for i in 1..<penalties.count { guard let left=penalties[i-1], penalties[i] == nil || left<penalties[i]! else { throw Failure.invalidInput("penalty order") } }
        var losses=Array(repeating:0.0,count:penalties.count)
        for gene in 0..<g {
            if gene%256==0 { try Task.checkCancellation() }
            let x=controls.map { $0[gene] }, d=zip(controls,treated).map { $1[gene]-$0[gene] }
            for holdout in 0..<n {
                let keep=(0..<n).filter { $0 != holdout }, xx=keep.map { x[$0] }, dd=keep.map { d[$0] }
                for k in penalties.indices {
                    let (xm,dm,b)=coefficients(xx,dd,penalties[k])
                    let predicted=max(0,x[holdout]+dm+b*(x[holdout]-xm))
                    let error=predicted-treated[holdout][gene];losses[k]+=error*error
                }
            }
        }
        losses=losses.map { $0/Double(n*g) }
        guard losses.allSatisfy({ $0.isFinite }) else { throw Failure.invalidInput("nonfinite validation loss") }
        var best=0
        for k in losses.indices where losses[k]<=losses[best] { best=k }
        var centers=[Double](),means=[Double](),slopes=[Double]()
        for gene in 0..<g {
            if gene%256==0 { try Task.checkCancellation() }
            let x=controls.map { $0[gene] },d=zip(controls,treated).map { $1[gene]-$0[gene] }
            let (xm,dm,b)=coefficients(x,d,penalties[best]);centers.append(xm);means.append(dm);slopes.append(b)
        }
        return Model(featureIDs:featureIDs,trainingDonorIDs:trainingDonorIDs,selectedPenalty:penalties[best],penalties:penalties,innerDonorMSE:losses,controlMean:centers,responseMean:means,responseSlope:slopes)
    }
    public static func predict(_ model: Model, featureIDs: [String], queryDonorID: String, control: [Double]) throws -> [Double] {
        let g=featureIDs.count
        guard g<=100_000,featureIDs==model.featureIDs,validIDs(featureIDs),validIDs(model.trainingDonorIDs),!queryDonorID.isEmpty,!model.trainingDonorIDs.contains(queryDonorID),control.count==g,
              control.allSatisfy({ $0.isFinite && $0>=0 && $0<=100 }),
              [model.controlMean,model.responseMean,model.responseSlope].allSatisfy({ $0.count==g && $0.allSatisfy(\.isFinite) }) else { throw Failure.invalidInput("query or model identity, dimensions or numeric values") }
        var values=[Double]()
        for i in 0..<g {
            let adjustment=model.responseSlope[i]*(control[i]-model.controlMean[i])
            values.append(max(0,control[i]+model.responseMean[i]+adjustment))
        }
        guard values.allSatisfy(\.isFinite) else { throw Failure.invalidInput("nonfinite prediction") }
        return values
    }
}
