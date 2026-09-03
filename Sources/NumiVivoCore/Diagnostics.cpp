#include "NumiVivoCore/Core.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace nvivo {

namespace {

void writeJsonString(std::ostringstream& stream, std::string_view value) {
    stream << '"' << json::escape(value) << '"';
}

} // namespace

void Diagnostics::add(Diagnostic diagnostic) {
    entries_.push_back(std::move(diagnostic));
}

void Diagnostics::note(std::string code,
                       std::string message,
                       std::string path,
                       std::string hint) {
    add({Severity::note, std::move(code), std::move(message), std::move(path), std::move(hint), {}});
}

void Diagnostics::warning(std::string code,
                          std::string message,
                          std::string path,
                          std::string hint) {
    add({Severity::warning, std::move(code), std::move(message), std::move(path), std::move(hint), {}});
}

void Diagnostics::error(std::string code,
                        std::string message,
                        std::string path,
                        std::string hint) {
    add({Severity::error, std::move(code), std::move(message), std::move(path), std::move(hint), {}});
}

void Diagnostics::fatal(std::string code,
                        std::string message,
                        std::string path,
                        std::string hint) {
    add({Severity::fatal, std::move(code), std::move(message), std::move(path), std::move(hint), {}});
}

bool Diagnostics::hasErrors() const noexcept {
    return std::any_of(entries_.begin(), entries_.end(), [](const Diagnostic& diagnostic) {
        return diagnostic.severity == Severity::error || diagnostic.severity == Severity::fatal;
    });
}

bool Diagnostics::hasFatal() const noexcept {
    return std::any_of(entries_.begin(), entries_.end(), [](const Diagnostic& diagnostic) {
        return diagnostic.severity == Severity::fatal;
    });
}

const std::vector<Diagnostic>& Diagnostics::entries() const noexcept {
    return entries_;
}

void Diagnostics::append(const Diagnostics& other) {
    entries_.insert(entries_.end(), other.entries_.begin(), other.entries_.end());
}

std::string Diagnostics::toJson() const {
    std::ostringstream stream;
    stream << "{\"diagnostics\":[";

    bool first = true;
    for (const auto& diagnostic : entries_) {
        if (!first) {
            stream << ',';
        }
        first = false;

        stream << "{\"severity\":";
        writeJsonString(stream, severityName(diagnostic.severity));
        stream << ",\"code\":";
        writeJsonString(stream, diagnostic.code);
        stream << ",\"message\":";
        writeJsonString(stream, diagnostic.message);
        stream << ",\"path\":";
        writeJsonString(stream, diagnostic.path);
        stream << ",\"hint\":";
        writeJsonString(stream, diagnostic.hint);
        stream << ",\"source\":{\"offset\":" << diagnostic.source.offset
               << ",\"length\":" << diagnostic.source.length
               << ",\"line\":" << diagnostic.source.line
               << ",\"column\":" << diagnostic.source.column << "}}";
    }

    stream << "],\"summary\":{\"count\":" << entries_.size()
           << ",\"hasErrors\":" << (hasErrors() ? "true" : "false")
           << ",\"hasFatal\":" << (hasFatal() ? "true" : "false")
           << "}}";
    return stream.str();
}

bool Dimension::isDimensionless() const noexcept {
    return length == 0 && mass == 0 && time == 0 && amount == 0 &&
           temperature == 0 && count == 0;
}

bool isValidIdentifier(std::string_view identifier) noexcept {
    if (identifier.empty() || identifier.size() > 128) {
        return false;
    }

    const auto first = static_cast<unsigned char>(identifier.front());
    if (!(std::isalpha(first) || identifier.front() == '_')) {
        return false;
    }

    return std::all_of(identifier.begin() + 1, identifier.end(), [](char character) {
        const auto value = static_cast<unsigned char>(character);
        return std::isalnum(value) || character == '_' || character == '-' || character == '.';
    });
}

std::string_view severityName(Severity severity) noexcept {
    switch (severity) {
        case Severity::note: return "note";
        case Severity::warning: return "warning";
        case Severity::error: return "error";
        case Severity::fatal: return "fatal";
    }
    return "unknown";
}

std::string_view fidelityName(FidelityLevel fidelity) noexcept {
    switch (fidelity) {
        case FidelityLevel::f0Logic: return "F0";
        case FidelityLevel::f1Deterministic: return "F1";
        case FidelityLevel::f2Stochastic: return "F2";
        case FidelityLevel::f3Spatial: return "F3";
        case FidelityLevel::f4Tissue: return "F4";
    }
    return "unknown";
}

std::string SafetyReport::toJson() const {
    std::ostringstream stream;
    stream << "{\"hasBlockingFinding\":" << (hasBlockingFinding ? "true" : "false")
           << ",\"monitoredOutputCount\":" << monitoredOutputCount
           << ",\"irreversibleActionCount\":" << irreversibleActionCount
           << ",\"findings\":[";

    bool first = true;
    for (const auto& finding : findings) {
        if (!first) {
            stream << ',';
        }
        first = false;

        stream << "{\"id\":";
        writeJsonString(stream, finding.id);
        stream << ",\"category\":" << static_cast<std::uint32_t>(finding.category)
               << ",\"severity\":";
        writeJsonString(stream, severityName(finding.severity));
        stream << ",\"subject\":";
        writeJsonString(stream, finding.subject);
        stream << ",\"message\":";
        writeJsonString(stream, finding.message);
        stream << ",\"requiredMitigation\":";
        writeJsonString(stream, finding.requiredMitigation);
        stream << '}';
    }

    stream << "]}";
    return stream.str();
}

} // namespace nvivo
