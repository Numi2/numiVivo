import Foundation

public struct VivoDimension: Hashable, Sendable, Codable {
    public var length: Int8
    public var mass: Int8
    public var time: Int8
    public var amount: Int8
    public var temperature: Int8
    public var count: Int8

    public init(
        length: Int8 = 0,
        mass: Int8 = 0,
        time: Int8 = 0,
        amount: Int8 = 0,
        temperature: Int8 = 0,
        count: Int8 = 0
    ) {
        self.length = length
        self.mass = mass
        self.time = time
        self.amount = amount
        self.temperature = temperature
        self.count = count
    }
}

public struct VivoUnitDefinition: Hashable, Sendable, Codable {
    public let symbol: String
    public let dimension: VivoDimension
    public let scaleToSI: Double
    public let offsetToSI: Double

    public init(symbol: String, dimension: VivoDimension, scaleToSI: Double, offsetToSI: Double = 0) {
        self.symbol = symbol
        self.dimension = dimension
        self.scaleToSI = scaleToSI
        self.offsetToSI = offsetToSI
    }
}

public struct VivoUnitSystem: Sendable {
    public static let standard = VivoUnitSystem()
    private let definitions: [String: VivoUnitDefinition]

    public init() {
        var definitions: [String: VivoUnitDefinition] = [:]
        func add(_ symbols: [String], dimension: VivoDimension, scale: Double, offset: Double = 0) {
            for symbol in symbols {
                definitions[symbol] = VivoUnitDefinition(
                    symbol: symbol,
                    dimension: dimension,
                    scaleToSI: scale,
                    offsetToSI: offset
                )
            }
        }

        let scalar = VivoDimension()
        let length = VivoDimension(length: 1)
        let area = VivoDimension(length: 2)
        let volume = VivoDimension(length: 3)
        let mass = VivoDimension(mass: 1)
        let time = VivoDimension(time: 1)
        let inverseTime = VivoDimension(time: -1)
        let amount = VivoDimension(amount: 1)
        let count = VivoDimension(count: 1)
        let concentration = VivoDimension(length: -3, amount: 1)
        let concentrationRate = VivoDimension(length: -3, time: -1, amount: 1)
        let diffusion = VivoDimension(length: 2, time: -1)
        let velocity = VivoDimension(length: 1, time: -1)
        let massRate = VivoDimension(mass: 1, time: -1)
        let countRate = VivoDimension(time: -1, count: 1)
        let secondOrderMolar = VivoDimension(length: 3, time: -1, amount: -1)
        let temperature = VivoDimension(temperature: 1)

        add(["1", "dimensionless", "normalized", "fraction"], dimension: scalar, scale: 1)
        add(["%"], dimension: scalar, scale: 0.01)

        add(["m"], dimension: length, scale: 1)
        add(["cm"], dimension: length, scale: 1e-2)
        add(["mm"], dimension: length, scale: 1e-3)
        add(["um", "µm"], dimension: length, scale: 1e-6)
        add(["nm"], dimension: length, scale: 1e-9)

        add(["m2", "m^2"], dimension: area, scale: 1)
        add(["cm2", "cm^2"], dimension: area, scale: 1e-4)
        add(["mm2", "mm^2"], dimension: area, scale: 1e-6)
        add(["um2", "um^2", "µm²"], dimension: area, scale: 1e-12)

        add(["m3", "m^3"], dimension: volume, scale: 1)
        add(["L"], dimension: volume, scale: 1e-3)
        add(["mL"], dimension: volume, scale: 1e-6)
        add(["uL", "µL"], dimension: volume, scale: 1e-9)
        add(["nL"], dimension: volume, scale: 1e-12)
        add(["pL"], dimension: volume, scale: 1e-15)
        add(["fL"], dimension: volume, scale: 1e-18)

        add(["kg"], dimension: mass, scale: 1)
        add(["g"], dimension: mass, scale: 1e-3)
        add(["mg"], dimension: mass, scale: 1e-6)
        add(["ug", "µg"], dimension: mass, scale: 1e-9)
        add(["ng"], dimension: mass, scale: 1e-12)
        add(["pg"], dimension: mass, scale: 1e-15)

        add(["s"], dimension: time, scale: 1)
        add(["ms"], dimension: time, scale: 1e-3)
        add(["us", "µs"], dimension: time, scale: 1e-6)
        add(["min"], dimension: time, scale: 60)
        add(["h", "hour"], dimension: time, scale: 3_600)
        add(["d", "day"], dimension: time, scale: 86_400)

        add(["Hz", "1/s", "s^-1"], dimension: inverseTime, scale: 1)
        add(["1/min", "min^-1"], dimension: inverseTime, scale: 1.0 / 60.0)
        add(["1/h", "h^-1"], dimension: inverseTime, scale: 1.0 / 3_600.0)
        add(["1/d", "d^-1"], dimension: inverseTime, scale: 1.0 / 86_400.0)

        add(["mol"], dimension: amount, scale: 1)
        add(["mmol"], dimension: amount, scale: 1e-3)
        add(["umol", "µmol"], dimension: amount, scale: 1e-6)
        add(["nmol"], dimension: amount, scale: 1e-9)
        add(["pmol"], dimension: amount, scale: 1e-12)

        add(["count", "molecule", "molecules", "copy", "copies"], dimension: count, scale: 1)

        add(["M"], dimension: concentration, scale: 1e3)
        add(["mM"], dimension: concentration, scale: 1)
        add(["uM", "µM"], dimension: concentration, scale: 1e-3)
        add(["nM"], dimension: concentration, scale: 1e-6)
        add(["pM"], dimension: concentration, scale: 1e-9)
        add(["fM"], dimension: concentration, scale: 1e-12)

        add(["M/s"], dimension: concentrationRate, scale: 1e3)
        add(["mM/s"], dimension: concentrationRate, scale: 1)
        add(["uM/s", "µM/s"], dimension: concentrationRate, scale: 1e-3)
        add(["nM/s"], dimension: concentrationRate, scale: 1e-6)
        add(["nM/min"], dimension: concentrationRate, scale: 1e-6 / 60)
        add(["nM/h"], dimension: concentrationRate, scale: 1e-6 / 3_600)

        add(["m2/s", "m^2/s"], dimension: diffusion, scale: 1)
        add(["cm2/s", "cm^2/s"], dimension: diffusion, scale: 1e-4)
        add(["mm2/s", "mm^2/s"], dimension: diffusion, scale: 1e-6)
        add(["um2/s", "um^2/s", "µm²/s"], dimension: diffusion, scale: 1e-12)

        add(["m/s"], dimension: velocity, scale: 1)
        add(["mm/s"], dimension: velocity, scale: 1e-3)
        add(["um/s", "µm/s"], dimension: velocity, scale: 1e-6)
        add(["um/min", "µm/min"], dimension: velocity, scale: 1e-6 / 60)

        add(["count/s", "molecule/s"], dimension: countRate, scale: 1)
        add(["count/min"], dimension: countRate, scale: 1.0 / 60.0)
        add(["ng/s"], dimension: massRate, scale: 1e-12)
        add(["ng/min"], dimension: massRate, scale: 1e-12 / 60)
        add(["ng/h", "ng/hour"], dimension: massRate, scale: 1e-12 / 3_600)

        add(["M^-1 s^-1", "1/(M*s)"], dimension: secondOrderMolar, scale: 1e-3)
        add(["mM^-1 s^-1"], dimension: secondOrderMolar, scale: 1)
        add(["uM^-1 s^-1", "µM^-1 s^-1"], dimension: secondOrderMolar, scale: 1e3)
        add(["nM^-1 s^-1", "1/(nM*s)"], dimension: secondOrderMolar, scale: 1e6)

        add(["K"], dimension: temperature, scale: 1)
        add(["degC", "°C"], dimension: temperature, scale: 1, offset: 273.15)

        self.definitions = definitions
    }

    public func definition(for symbol: String) -> VivoUnitDefinition? {
        definitions[symbol]
    }

    public func areCompatible(_ lhs: String, _ rhs: String) -> Bool {
        guard let left = definitions[lhs], let right = definitions[rhs] else { return false }
        return left.dimension == right.dimension
    }

    public func convert(_ value: Double, from source: String, to destination: String) throws -> Double {
        guard value.isFinite else {
            throw VivoArtifactValidationError.invalid("unit conversion input must be finite")
        }
        guard let sourceDefinition = definitions[source] else {
            throw VivoArtifactValidationError.unresolved("unknown unit \(source)")
        }
        guard let destinationDefinition = definitions[destination] else {
            throw VivoArtifactValidationError.unresolved("unknown unit \(destination)")
        }
        guard sourceDefinition.dimension == destinationDefinition.dimension else {
            throw VivoArtifactValidationError.incompatible("units \(source) and \(destination) have different dimensions")
        }
        let si = (value + sourceDefinition.offsetToSI) * sourceDefinition.scaleToSI
        let converted = si / destinationDefinition.scaleToSI - destinationDefinition.offsetToSI
        guard converted.isFinite else {
            throw VivoArtifactValidationError.invalid("unit conversion produced a non-finite value")
        }
        return converted
    }

    public func convert(_ quantity: VivoQuantity, to destination: String) throws -> Double {
        try convert(quantity.value, from: quantity.unit, to: destination)
    }
}
