#include "NumiVivoCore/ScalableReactionNetworkAnalysis.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace nvivo::network {

namespace {

struct RREFResult {
    std::vector<std::vector<double>> matrix;
    std::vector<std::uint32_t> pivotColumns;
    std::uint32_t rank = 0;
    double minimumAcceptedPivot = std::numeric_limits<double>::infinity();
    double maximumInputMagnitude = 0.0;
};

RREFResult rref(std::vector<std::vector<double>> matrix,
                double tolerance) {
    RREFResult result;
    result.matrix = std::move(matrix);
    const std::size_t rows = result.matrix.size();
    const std::size_t columns = rows == 0 ? 0 : result.matrix.front().size();

    for (const auto& row : result.matrix) {
        for (const double value : row) {
            result.maximumInputMagnitude = std::max(
                result.maximumInputMagnitude,
                std::abs(value)
            );
        }
    }
    const double absoluteTolerance = tolerance *
        std::max(1.0, result.maximumInputMagnitude);

    std::size_t pivotRow = 0;
    for (std::size_t column = 0;
         column < columns && pivotRow < rows;
         ++column) {
        std::size_t selected = pivotRow;
        double selectedMagnitude = 0.0;
        for (std::size_t row = pivotRow; row < rows; ++row) {
            const double magnitude = std::abs(result.matrix[row][column]);
            if (magnitude > selectedMagnitude) {
                selected = row;
                selectedMagnitude = magnitude;
            }
        }
        if (selectedMagnitude <= absoluteTolerance) continue;
        if (selected != pivotRow) {
            std::swap(result.matrix[selected], result.matrix[pivotRow]);
        }

        const double pivot = result.matrix[pivotRow][column];
        result.minimumAcceptedPivot = std::min(
            result.minimumAcceptedPivot,
            std::abs(pivot)
        );
        for (std::size_t entry = column; entry < columns; ++entry) {
            result.matrix[pivotRow][entry] /= pivot;
            if (std::abs(result.matrix[pivotRow][entry]) <= absoluteTolerance) {
                result.matrix[pivotRow][entry] = 0.0;
            }
        }

        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivotRow) continue;
            const double factor = result.matrix[row][column];
            if (std::abs(factor) <= absoluteTolerance) {
                result.matrix[row][column] = 0.0;
                continue;
            }
            for (std::size_t entry = column; entry < columns; ++entry) {
                result.matrix[row][entry] -=
                    factor * result.matrix[pivotRow][entry];
                if (std::abs(result.matrix[row][entry]) <= absoluteTolerance) {
                    result.matrix[row][entry] = 0.0;
                }
            }
        }
        result.pivotColumns.push_back(static_cast<std::uint32_t>(column));
        ++pivotRow;
    }
    result.rank = static_cast<std::uint32_t>(pivotRow);
    return result;
}

bool componentFits(const Component& component,
                   const Options& options) {
    const std::uint64_t species = component.speciesIndices.size();
    const std::uint64_t reactions = component.reactionIndices.size();
    const auto product = species * reactions;
    return species <= options.maximumDenseSpecies &&
           reactions <= options.maximumDenseReactions &&
           product <= options.maximumDenseElements;
}

bool buildStoichiometry(
    const ProgramIR& ir,
    const Component& component,
    std::vector<std::vector<double>>& stoichiometry) {
    const std::size_t speciesCount = component.speciesIndices.size();
    const std::size_t reactionCount = component.reactionIndices.size();
    stoichiometry.assign(
        speciesCount,
        std::vector<double>(reactionCount, 0.0)
    );

    std::unordered_map<std::uint32_t, std::uint32_t> localSpecies;
    localSpecies.reserve(speciesCount);
    for (std::uint32_t local = 0; local < speciesCount; ++local) {
        localSpecies.emplace(component.speciesIndices[local], local);
    }

    for (std::uint32_t localReaction = 0;
         localReaction < reactionCount;
         ++localReaction) {
        const std::uint32_t reactionIndex =
            component.reactionIndices[localReaction];
        if (reactionIndex >= ir.reactions.size()) return false;
        const ReactionRecord& reaction = ir.reactions[reactionIndex];
        const std::uint64_t reactantEnd =
            static_cast<std::uint64_t>(reaction.reactantOffset) +
            reaction.reactantCount;
        const std::uint64_t productEnd =
            static_cast<std::uint64_t>(reaction.productOffset) +
            reaction.productCount;
        if (reactantEnd > ir.stoichiometry.size() ||
            productEnd > ir.stoichiometry.size()) {
            return false;
        }

        for (std::uint32_t termIndex = 0;
             termIndex < reaction.reactantCount;
             ++termIndex) {
            const auto& term =
                ir.stoichiometry[reaction.reactantOffset + termIndex];
            const auto found = localSpecies.find(term.speciesIndex);
            if (found == localSpecies.end() || term.coefficient <= 0) {
                return false;
            }
            stoichiometry[found->second][localReaction] -= term.coefficient;
        }
        for (std::uint32_t termIndex = 0;
             termIndex < reaction.productCount;
             ++termIndex) {
            const auto& term =
                ir.stoichiometry[reaction.productOffset + termIndex];
            const auto found = localSpecies.find(term.speciesIndex);
            if (found == localSpecies.end() || term.coefficient <= 0) {
                return false;
            }
            stoichiometry[found->second][localReaction] += term.coefficient;
        }
    }
    return true;
}

std::vector<std::vector<double>> transpose(
    const std::vector<std::vector<double>>& matrix,
    std::size_t columns) {
    const std::size_t rows = matrix.size();
    std::vector<std::vector<double>> result(
        columns,
        std::vector<double>(rows, 0.0)
    );
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t column = 0; column < columns; ++column) {
            result[column][row] = matrix[row][column];
        }
    }
    return result;
}

std::vector<ConservationLaw> inferConservation(
    const RREFResult& decomposition,
    const std::vector<std::vector<double>>& stoichiometry,
    const Component& component,
    double tolerance,
    std::uint32_t maximumLaws) {
    const std::size_t speciesCount = component.speciesIndices.size();
    const std::size_t reactionCount = component.reactionIndices.size();
    std::set<std::uint32_t> pivotSet(
        decomposition.pivotColumns.begin(),
        decomposition.pivotColumns.end()
    );
    std::vector<ConservationLaw> laws;

    for (std::uint32_t freeColumn = 0;
         freeColumn < speciesCount && laws.size() < maximumLaws;
         ++freeColumn) {
        if (pivotSet.contains(freeColumn)) continue;

        std::vector<double> vector(speciesCount, 0.0);
        vector[freeColumn] = 1.0;
        for (std::size_t row = 0;
             row < decomposition.pivotColumns.size();
             ++row) {
            const std::uint32_t pivot = decomposition.pivotColumns[row];
            vector[pivot] = -decomposition.matrix[row][freeColumn];
        }

        double maximumMagnitude = 0.0;
        for (const double value : vector) {
            maximumMagnitude = std::max(maximumMagnitude, std::abs(value));
        }
        if (!(maximumMagnitude > 0.0) || !std::isfinite(maximumMagnitude)) {
            continue;
        }
        for (double& value : vector) value /= maximumMagnitude;
        const auto first = std::find_if(
            vector.begin(),
            vector.end(),
            [tolerance](double value) {
                return std::abs(value) > tolerance;
            }
        );
        if (first != vector.end() && *first < 0.0) {
            for (double& value : vector) value = -value;
        }

        ConservationLaw law;
        for (std::uint32_t localSpecies = 0;
             localSpecies < speciesCount;
             ++localSpecies) {
            if (std::abs(vector[localSpecies]) > tolerance) {
                law.terms.push_back({
                    component.speciesIndices[localSpecies],
                    vector[localSpecies]
                });
            }
        }
        if (law.terms.size() < 2) continue;

        for (std::size_t localReaction = 0;
             localReaction < reactionCount;
             ++localReaction) {
            double residual = 0.0;
            for (std::size_t localSpecies = 0;
                 localSpecies < speciesCount;
                 ++localSpecies) {
                residual += vector[localSpecies] *
                    stoichiometry[localSpecies][localReaction];
            }
            law.maximumResidual = std::max(
                law.maximumResidual,
                std::abs(residual)
            );
        }
        if (law.maximumResidual <= tolerance * 10.0) {
            laws.push_back(std::move(law));
        }
    }
    return laws;
}

bool nearTolerance(const RREFResult& decomposition,
                   double rankTolerance) {
    return decomposition.rank > 0 &&
           decomposition.minimumAcceptedPivot <=
               rankTolerance *
               std::max(1.0, decomposition.maximumInputMagnitude) *
               100.0;
}

} // namespace

Result ScalableAnalyzer::analyze(const ProgramIR& ir,
                                 const Options& options) const {
    Result result = Analyzer().analyze(ir, options);
    if (result.denseAnalysisComplete || result.diagnostics.hasErrors()) {
        return result;
    }

    std::uint64_t totalRank = 0;
    std::vector<ConservationLaw> laws;
    bool complete = true;
    bool uncertain = false;
    std::uint32_t oversizedComponent = 0;

    for (std::uint32_t componentIndex = 0;
         componentIndex < result.components.size();
         ++componentIndex) {
        const Component& component = result.components[componentIndex];
        if (!componentFits(component, options)) {
            complete = false;
            oversizedComponent = componentIndex;
            break;
        }

        std::vector<std::vector<double>> stoichiometry;
        if (!buildStoichiometry(ir, component, stoichiometry)) {
            complete = false;
            oversizedComponent = componentIndex;
            break;
        }

        const auto transposed = transpose(
            stoichiometry,
            component.reactionIndices.size()
        );
        const RREFResult decomposition = rref(
            transposed,
            options.rankTolerance
        );
        totalRank += decomposition.rank;
        uncertain = uncertain || nearTolerance(
            decomposition,
            options.rankTolerance
        );

        if (options.inferConservationLaws &&
            laws.size() < options.maximumConservationLaws) {
            const std::uint32_t remaining =
                options.maximumConservationLaws -
                static_cast<std::uint32_t>(laws.size());
            auto componentLaws = inferConservation(
                decomposition,
                stoichiometry,
                component,
                options.conservationTolerance,
                remaining
            );
            laws.insert(
                laws.end(),
                std::make_move_iterator(componentLaws.begin()),
                std::make_move_iterator(componentLaws.end())
            );
        }
    }

    if (!complete || totalRank > std::numeric_limits<std::uint32_t>::max()) {
        result.issues.push_back({
            IssueCategory::denseAnalysisSkipped,
            Severity::note,
            oversizedComponent,
            "Componentwise rank/conservation fallback could not complete because at least one connected component exceeded the bounded dense-analysis budget or had invalid topology."
        });
        return result;
    }

    result.stoichiometricRank = static_cast<std::uint32_t>(totalRank);
    result.denseAnalysisComplete = true;
    if (options.inferConservationLaws) {
        std::sort(
            laws.begin(),
            laws.end(),
            [](const ConservationLaw& left, const ConservationLaw& right) {
                const auto leftFirst = left.terms.empty()
                    ? std::numeric_limits<std::uint32_t>::max()
                    : left.terms.front().speciesIndex;
                const auto rightFirst = right.terms.empty()
                    ? std::numeric_limits<std::uint32_t>::max()
                    : right.terms.front().speciesIndex;
                if (leftFirst != rightFirst) return leftFirst < rightFirst;
                return left.terms.size() < right.terms.size();
            }
        );
        result.conservationLaws = std::move(laws);
    }

    result.issues.erase(
        std::remove_if(
            result.issues.begin(),
            result.issues.end(),
            [](const Issue& issue) {
                return issue.category == IssueCategory::denseAnalysisSkipped;
            }
        ),
        result.issues.end()
    );
    result.issues.push_back({
        IssueCategory::denseAnalysisSkipped,
        Severity::note,
        0,
        "Global dense rank/conservation analysis exceeded its matrix budget; exact analysis completed componentwise across disconnected reaction-network components."
    });

    if (uncertain) {
        const bool alreadyPresent = std::any_of(
            result.issues.begin(),
            result.issues.end(),
            [](const Issue& issue) {
                return issue.category == IssueCategory::numericalRankUncertain;
            }
        );
        if (!alreadyPresent) {
            result.issues.push_back({
                IssueCategory::numericalRankUncertain,
                Severity::warning,
                0,
                "At least one componentwise stoichiometric-rank pivot is close to the configured tolerance."
            });
        }
    }

    if (options.computeComplexDeficiency) {
        result.deficiency =
            static_cast<std::int64_t>(result.complexCount) -
            static_cast<std::int64_t>(result.linkageClassCount) -
            static_cast<std::int64_t>(result.stoichiometricRank);
        const bool hasPositiveDeficiency = std::any_of(
            result.issues.begin(),
            result.issues.end(),
            [](const Issue& issue) {
                return issue.category == IssueCategory::positiveDeficiency;
            }
        );
        if (result.deficiency < 0) {
            result.issues.push_back({
                IssueCategory::numericalRankUncertain,
                Severity::warning,
                0,
                "Componentwise rank produced a negative global deficiency; topology or numerical rank requires review."
            });
        } else if (result.deficiency > 0 && !hasPositiveDeficiency) {
            result.issues.push_back({
                IssueCategory::positiveDeficiency,
                Severity::note,
                0,
                "Reaction network has positive structural deficiency " +
                    std::to_string(result.deficiency) +
                    " after componentwise rank analysis."
            });
        }
    }

    return result;
}

} // namespace nvivo::network
