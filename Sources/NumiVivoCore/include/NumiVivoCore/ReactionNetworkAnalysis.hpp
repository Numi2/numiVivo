#pragma once

#include "NumiVivoCore/Core.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace nvivo::network {

enum class IssueCategory : std::uint8_t {
    isolatedSpecies = 0,
    sourceOnlySpecies = 1,
    sinkOnlySpecies = 2,
    constitutiveSource = 3,
    terminalSink = 4,
    declaredConservationViolation = 5,
    duplicateReactionVector = 6,
    emptyReaction = 7,
    nonWeaklyReversible = 8,
    positiveDeficiency = 9,
    denseAnalysisSkipped = 10,
    numericalRankUncertain = 11,
    invalidTopology = 12
};

struct Issue {
    IssueCategory category = IssueCategory::invalidTopology;
    Severity severity = Severity::warning;
    std::uint32_t subjectIndex = 0;
    std::string message;
};

struct ConservationTerm {
    std::uint32_t speciesIndex = 0;
    double coefficient = 0.0;
};

struct ConservationLaw {
    std::vector<ConservationTerm> terms;
    double maximumResidual = 0.0;
};

struct Component {
    std::vector<std::uint32_t> speciesIndices;
    std::vector<std::uint32_t> reactionIndices;
};

struct ComplexRecord {
    std::vector<StoichiometryRecord> terms;
    std::vector<std::uint32_t> outgoingComplexes;
    std::vector<std::uint32_t> incomingComplexes;
};

struct Options {
    double rankTolerance = 1.0e-10;
    double conservationTolerance = 1.0e-9;
    std::uint32_t maximumDenseSpecies = 2'048;
    std::uint32_t maximumDenseReactions = 2'048;
    std::uint64_t maximumDenseElements = 16'777'216;
    std::uint32_t maximumConservationLaws = 256;
    bool inferConservationLaws = true;
    bool computeComplexDeficiency = true;
};

struct Result {
    std::uint32_t speciesCount = 0;
    std::uint32_t reactionCount = 0;
    std::uint32_t complexCount = 0;
    std::uint32_t linkageClassCount = 0;
    std::uint32_t strongLinkageClassCount = 0;
    std::uint32_t stoichiometricRank = 0;
    std::int64_t deficiency = 0;
    bool weaklyReversible = false;
    bool denseAnalysisComplete = false;
    std::vector<Component> components;
    std::vector<ComplexRecord> complexes;
    std::vector<ConservationLaw> conservationLaws;
    std::vector<Issue> issues;
    Diagnostics diagnostics;

    [[nodiscard]] std::string toJson() const;
};

class Analyzer {
public:
    [[nodiscard]] Result analyze(const ProgramIR& ir,
                                 const Options& options = {}) const;
};

[[nodiscard]] std::string_view issueCategoryName(IssueCategory category) noexcept;

} // namespace nvivo::network
