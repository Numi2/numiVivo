#include "NumiVivoCore/Core.hpp"

#include <cerrno>
#include <cmath>
#include <cstdlib>

namespace nvivo {

namespace {

std::string trim(std::string_view value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) {
        return {};
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return std::string(value.substr(first, last - first + 1));
}

Dimension makeDimension(std::int8_t length,
                        std::int8_t mass,
                        std::int8_t time,
                        std::int8_t amount,
                        std::int8_t temperature,
                        std::int8_t count) {
    return {length, mass, time, amount, temperature, count};
}

} // namespace

UnitRegistry::UnitRegistry() {
    const Dimension scalar{};
    const Dimension length = makeDimension(1, 0, 0, 0, 0, 0);
    const Dimension area = makeDimension(2, 0, 0, 0, 0, 0);
    const Dimension volume = makeDimension(3, 0, 0, 0, 0, 0);
    const Dimension mass = makeDimension(0, 1, 0, 0, 0, 0);
    const Dimension time = makeDimension(0, 0, 1, 0, 0, 0);
    const Dimension inverseTime = makeDimension(0, 0, -1, 0, 0, 0);
    const Dimension amount = makeDimension(0, 0, 0, 1, 0, 0);
    const Dimension discreteCount = makeDimension(0, 0, 0, 0, 0, 1);
    const Dimension concentration = makeDimension(-3, 0, 0, 1, 0, 0);
    const Dimension concentrationRate = makeDimension(-3, 0, -1, 1, 0, 0);
    const Dimension diffusion = makeDimension(2, 0, -1, 0, 0, 0);
    const Dimension countRate = makeDimension(0, 0, -1, 0, 0, 1);
    const Dimension massRate = makeDimension(0, 1, -1, 0, 0, 0);
    const Dimension secondOrderMolar = makeDimension(3, 0, -1, -1, 0, 0);

    const auto add = [this](std::string symbol,
                            Dimension valueDimension,
                            double scale,
                            bool contextDependent = false) {
        UnitDefinition definition{symbol, valueDimension, scale, 0.0, contextDependent};
        units_.emplace(std::move(symbol), std::move(definition));
    };

    add("1", scalar, 1.0);
    add("dimensionless", scalar, 1.0);
    add("normalized", scalar, 1.0);
    add("fraction", scalar, 1.0);
    add("%", scalar, 0.01);

    add("m", length, 1.0);
    add("cm", length, 1.0e-2);
    add("mm", length, 1.0e-3);
    add("um", length, 1.0e-6);
    add("µm", length, 1.0e-6);
    add("nm", length, 1.0e-9);
    add("m2", area, 1.0);
    add("um2", area, 1.0e-12);
    add("µm²", area, 1.0e-12);

    add("m3", volume, 1.0);
    add("L", volume, 1.0e-3);
    add("mL", volume, 1.0e-6);
    add("uL", volume, 1.0e-9);
    add("µL", volume, 1.0e-9);
    add("nL", volume, 1.0e-12);
    add("pL", volume, 1.0e-15);
    add("fL", volume, 1.0e-18);

    add("kg", mass, 1.0);
    add("g", mass, 1.0e-3);
    add("mg", mass, 1.0e-6);
    add("ug", mass, 1.0e-9);
    add("µg", mass, 1.0e-9);
    add("ng", mass, 1.0e-12);
    add("pg", mass, 1.0e-15);

    add("s", time, 1.0);
    add("ms", time, 1.0e-3);
    add("us", time, 1.0e-6);
    add("µs", time, 1.0e-6);
    add("min", time, 60.0);
    add("h", time, 3'600.0);
    add("hour", time, 3'600.0);
    add("d", time, 86'400.0);
    add("day", time, 86'400.0);

    add("Hz", inverseTime, 1.0);
    add("1/s", inverseTime, 1.0);
    add("s^-1", inverseTime, 1.0);
    add("1/min", inverseTime, 1.0 / 60.0);
    add("min^-1", inverseTime, 1.0 / 60.0);
    add("1/h", inverseTime, 1.0 / 3'600.0);
    add("h^-1", inverseTime, 1.0 / 3'600.0);

    add("mol", amount, 1.0);
    add("mmol", amount, 1.0e-3);
    add("umol", amount, 1.0e-6);
    add("µmol", amount, 1.0e-6);
    add("nmol", amount, 1.0e-9);
    add("pmol", amount, 1.0e-12);

    add("count", discreteCount, 1.0);
    add("molecule", discreteCount, 1.0);
    add("copy", discreteCount, 1.0);

    // Concentration scales use mol / cubic metre internally.
    add("M", concentration, 1.0e3);
    add("mM", concentration, 1.0);
    add("uM", concentration, 1.0e-3);
    add("µM", concentration, 1.0e-3);
    add("nM", concentration, 1.0e-6);
    add("pM", concentration, 1.0e-9);
    add("fM", concentration, 1.0e-12);

    add("M/s", concentrationRate, 1.0e3);
    add("mM/s", concentrationRate, 1.0);
    add("uM/s", concentrationRate, 1.0e-3);
    add("µM/s", concentrationRate, 1.0e-3);
    add("nM/s", concentrationRate, 1.0e-6);
    add("nM/min", concentrationRate, 1.0e-6 / 60.0);

    add("m2/s", diffusion, 1.0);
    add("um2/s", diffusion, 1.0e-12);
    add("µm²/s", diffusion, 1.0e-12);

    add("count/s", countRate, 1.0);
    add("count/min", countRate, 1.0 / 60.0);
    add("ng/s", massRate, 1.0e-12);
    add("ng/min", massRate, 1.0e-12 / 60.0);
    add("ng/hour", massRate, 1.0e-12 / 3'600.0);

    add("M^-1 s^-1", secondOrderMolar, 1.0e-3);
    add("mM^-1 s^-1", secondOrderMolar, 1.0);
    add("uM^-1 s^-1", secondOrderMolar, 1.0e3);
    add("nM^-1 s^-1", secondOrderMolar, 1.0e6);
}

const UnitDefinition* UnitRegistry::find(std::string_view symbol) const noexcept {
    const auto iterator = units_.find(symbol);
    return iterator == units_.end() ? nullptr : &iterator->second;
}

bool UnitRegistry::compatible(std::string_view lhs, std::string_view rhs) const noexcept {
    const auto* left = find(lhs);
    const auto* right = find(rhs);
    return left != nullptr && right != nullptr && left->dimension == right->dimension;
}

std::optional<double> UnitRegistry::convert(double value,
                                            std::string_view from,
                                            std::string_view to) const noexcept {
    const auto* source = find(from);
    const auto* destination = find(to);
    if (source == nullptr || destination == nullptr || source->dimension != destination->dimension ||
        source->contextDependent != destination->contextDependent) {
        return std::nullopt;
    }

    const double si = value * source->scaleToSI;
    const double converted = si / destination->scaleToSI;
    return std::isfinite(converted) ? std::optional<double>(converted) : std::nullopt;
}

std::optional<double> UnitRegistry::parseDurationSeconds(std::string_view text) const noexcept {
    const std::string input = trim(text);
    if (input.empty()) {
        return std::nullopt;
    }

    errno = 0;
    char* end = nullptr;
    const double value = std::strtod(input.c_str(), &end);
    if (errno == ERANGE || end == input.c_str() || !std::isfinite(value) || value < 0.0) {
        return std::nullopt;
    }

    const std::string unitText = trim(std::string_view(end, input.c_str() + input.size() - end));
    const auto* unit = find(unitText);
    if (unit == nullptr || unit->dimension != makeDimension(0, 0, 1, 0, 0, 0)) {
        return std::nullopt;
    }

    const double seconds = value * unit->scaleToSI;
    return std::isfinite(seconds) ? std::optional<double>(seconds) : std::nullopt;
}

} // namespace nvivo
