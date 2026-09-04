#include "NumiVivoCore/ReactionNetworkAnalysis.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace nvivo::network {

namespace {

class UnionFind {
public:
    explicit UnionFind(std::size_t count) : parent_(count), rank_(count, 0) {
        std::iota(parent_.begin(), parent_.end(), std::size_t{0});
    }

    std::size_t find(std::size_t value) {
        if (parent_[value] != value) parent_[value] = find(parent_[value]);
        return parent_[value];
    }

    void unite(std::size_t left, std::size_t right) {
        left = find(left);
        right = find(right);
        if (left == right) return;
        if (rank_[left] < rank_[right]) {
            parent_[left] = right;
        } else if (rank_[left] > rank_[right]) {
            parent_[right] = left;
        } else {
            parent_[right] = left;
            ++rank_[left];
        }
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<std::uint8_t> rank_;
};

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
            result.maximumInputMagnitude = std::max(result.maximumInputMagnitude, std::abs(value));
        }
    }
    const double scale = std::max(1.0, result.maximumInputMagnitude);
    const double absoluteTolerance = tolerance * scale;

    std::size_t pivotRow = 0;
    for (std::size_t column = 0; column < columns && pivotRow < rows; ++column) {
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
        if (selected != pivotRow) std::swap(result.matrix[selected], result.matrix[pivotRow]);

        const double pivot = result.matrix[pivotRow][column];
        result.minimumAcceptedPivot = std::min(result.minimumAcceptedPivot, std::abs(pivot));
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
                result.matrix[row][entry] -= factor * result.matrix[pivotRow][entry];
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

std::vector<ConservationLaw> inferConservation(
    const RREFResult& decomposition,
    const std::vector<std::vector<double>>& stoichiometry,
    double tolerance,
    std::uint32_t maximumLaws) {
    const std::size_t speciesCount = decomposition.matrix.empty()
        ? (stoichiometry.empty() ? 0 : stoichiometry.size())
        : decomposition.matrix.front().size();
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
        for (std::size_t row = 0; row < decomposition.pivotColumns.size(); ++row) {
            const std::uint32_t pivot = decomposition.pivotColumns[row];
            vector[pivot] = -decomposition.matrix[row][freeColumn];
        }

        double maximum = 0.0;
        for (const double value : vector) maximum = std::max(maximum, std::abs(value));
        if (!(maximum > 0.0) || !std::isfinite(maximum)) continue;
        for (double& value : vector) value /= maximum;
        const auto first = std::find_if(vector.begin(), vector.end(), [tolerance](double value) {
            return std::abs(value) > tolerance;
        });
        if (first != vector.end() && *first < 0.0) {
            for (double& value : vector) value = -value;
        }

        ConservationLaw law;
        for (std::uint32_t species = 0; species < vector.size(); ++species) {
            if (std::abs(vector[species]) > tolerance) {
                law.terms.push_back({species, vector[species]});
            }
        }
        if (law.terms.size() < 2) continue;

        const std::size_t reactionCount = stoichiometry.empty() ? 0 : stoichiometry.front().size();
        for (std::size_t reaction = 0; reaction < reactionCount; ++reaction) {
            double residual = 0.0;
            for (std::size_t species = 0; species < speciesCount; ++species) {
                residual += vector[species] * stoichiometry[species][reaction];
            }
            law.maximumResidual = std::max(law.maximumResidual, std::abs(residual));
        }
        if (law.maximumResidual <= tolerance * 10.0) laws.push_back(std::move(law));
    }

    std::sort(laws.begin(), laws.end(), [](const ConservationLaw& left, const ConservationLaw& right) {
        const std::uint32_t leftFirst = left.terms.empty() ? std::numeric_limits<std::uint32_t>::max() : left.terms.front().speciesIndex;
        const std::uint32_t rightFirst = right.terms.empty() ? std::numeric_limits<std::uint32_t>::max() : right.terms.front().speciesIndex;
        if (leftFirst != rightFirst) return leftFirst < rightFirst;
        return left.terms.size() < right.terms.size();
    });
    return laws;
}

std::vector<StoichiometryRecord> canonicalComplex(
    std::span<const StoichiometryRecord> terms) {
    std::map<std::uint32_t, std::int32_t> combined;
    for (const auto& term : terms) {
        combined[term.speciesIndex] += std::abs(static_cast<std::int32_t>(term.coefficient));
    }
    std::vector<StoichiometryRecord> result;
    result.reserve(combined.size());
    for (const auto& [species, coefficient] : combined) {
        if (coefficient <= 0 || coefficient > std::numeric_limits<std::int16_t>::max()) continue;
        result.push_back({species, static_cast<std::int16_t>(coefficient), 0});
    }
    return result;
}

std::string complexKey(const std::vector<StoichiometryRecord>& terms) {
    std::ostringstream stream;
    for (const auto& term : terms) {
        stream << term.speciesIndex << ':' << term.coefficient << ';';
    }
    return stream.str();
}

std::string reactionVectorKey(const std::map<std::uint32_t, std::int32_t>& net) {
    std::ostringstream stream;
    for (const auto& [species, coefficient] : net) {
        if (coefficient != 0) stream << species << ':' << coefficient << ';';
    }
    return stream.str();
}

void tarjanVisit(
    std::uint32_t node,
    const std::vector<ComplexRecord>& complexes,
    std::vector<std::int32_t>& indices,
    std::vector<std::int32_t>& lowLink,
    std::vector<bool>& onStack,
    std::vector<std::uint32_t>& stack,
    std::int32_t& nextIndex,
    std::vector<std::uint32_t>& componentByNode,
    std::uint32_t& componentCount) {
    indices[node] = nextIndex;
    lowLink[node] = nextIndex;
    ++nextIndex;
    stack.push_back(node);
    onStack[node] = true;

    for (const std::uint32_t destination : complexes[node].outgoingComplexes) {
        if (indices[destination] < 0) {
            tarjanVisit(
                destination,
                complexes,
                indices,
                lowLink,
                onStack,
                stack,
                nextIndex,
                componentByNode,
                componentCount
            );
            lowLink[node] = std::min(lowLink[node], lowLink[destination]);
        } else if (onStack[destination]) {
            lowLink[node] = std::min(lowLink[node], indices[destination]);
        }
    }

    if (lowLink[node] == indices[node]) {
        while (!stack.empty()) {
            const std::uint32_t member = stack.back();
            stack.pop_back();
            onStack[member] = false;
            componentByNode[member] = componentCount;
            if (member == node) break;
        }
        ++componentCount;
    }
}

void addIssue(Result& result,
              IssueCategory category,
              Severity severity,
              std::uint32_t subject,
              std::string message) {
    result.issues.push_back({category, severity, subject, std::move(message)});
}

void writeJsonString(std::ostringstream& stream, std::string_view value) {
    stream << '"' << json::escape(value) << '"';
}

} // namespace

Result Analyzer::analyze(const ProgramIR& ir,
                         const Options& options) const {
    Result result;
    result.speciesCount = static_cast<std::uint32_t>(ir.species.size());
    result.reactionCount = static_cast<std::uint32_t>(ir.reactions.size());

    if (!(options.rankTolerance > 0.0) || !std::isfinite(options.rankTolerance) ||
        !(options.conservationTolerance > 0.0) || !std::isfinite(options.conservationTolerance) ||
        options.maximumDenseSpecies == 0 || options.maximumDenseReactions == 0 ||
        options.maximumDenseElements == 0) {
        result.diagnostics.error("NVN001", "Reaction-network analysis options are invalid.");
        return result;
    }
    if (ir.species.size() > std::numeric_limits<std::uint32_t>::max() ||
        ir.reactions.size() > std::numeric_limits<std::uint32_t>::max()) {
        result.diagnostics.error("NVN002", "Reaction network exceeds the UInt32 analysis ABI.");
        return result;
    }

    const std::size_t speciesCount = ir.species.size();
    const std::size_t reactionCount = ir.reactions.size();
    std::vector<bool> consumed(speciesCount, false);
    std::vector<bool> produced(speciesCount, false);
    UnionFind networkComponents(speciesCount + reactionCount);
    std::map<std::string, std::uint32_t, std::less<>> complexIDs;
    std::map<std::string, std::uint32_t, std::less<>> reactionVectors;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> complexEdges;

    auto internComplex = [&](std::vector<StoichiometryRecord> terms) -> std::uint32_t {
        const std::string key = complexKey(terms);
        const auto found = complexIDs.find(key);
        if (found != complexIDs.end()) return found->second;
        const auto index = static_cast<std::uint32_t>(result.complexes.size());
        complexIDs.emplace(key, index);
        result.complexes.push_back({std::move(terms), {}, {}});
        return index;
    };

    for (std::uint32_t reactionIndex = 0; reactionIndex < reactionCount; ++reactionIndex) {
        const ReactionRecord& reaction = ir.reactions[reactionIndex];
        const std::uint64_t reactantEnd = static_cast<std::uint64_t>(reaction.reactantOffset) + reaction.reactantCount;
        const std::uint64_t productEnd = static_cast<std::uint64_t>(reaction.productOffset) + reaction.productCount;
        if (reactantEnd > ir.stoichiometry.size() || productEnd > ir.stoichiometry.size()) {
            addIssue(
                result,
                IssueCategory::invalidTopology,
                Severity::error,
                reactionIndex,
                "Reaction stoichiometry range is outside the compiled table."
            );
            continue;
        }

        const auto reactants = std::span<const StoichiometryRecord>(ir.stoichiometry).subspan(
            reaction.reactantOffset,
            reaction.reactantCount
        );
        const auto products = std::span<const StoichiometryRecord>(ir.stoichiometry).subspan(
            reaction.productOffset,
            reaction.productCount
        );
        std::map<std::uint32_t, std::int32_t> net;
        bool topologyValid = true;
        for (const auto& term : reactants) {
            if (term.speciesIndex >= speciesCount || term.coefficient <= 0) {
                topologyValid = false;
                continue;
            }
            consumed[term.speciesIndex] = true;
            net[term.speciesIndex] -= term.coefficient;
            networkComponents.unite(term.speciesIndex, speciesCount + reactionIndex);
        }
        for (const auto& term : products) {
            if (term.speciesIndex >= speciesCount || term.coefficient <= 0) {
                topologyValid = false;
                continue;
            }
            produced[term.speciesIndex] = true;
            net[term.speciesIndex] += term.coefficient;
            networkComponents.unite(term.speciesIndex, speciesCount + reactionIndex);
        }
        if (!topologyValid) {
            addIssue(
                result,
                IssueCategory::invalidTopology,
                Severity::error,
                reactionIndex,
                "Reaction contains an invalid species index or coefficient."
            );
            continue;
        }
        if (reactants.empty() && products.empty()) {
            addIssue(
                result,
                IssueCategory::emptyReaction,
                Severity::error,
                reactionIndex,
                "Reaction has no reactants and no products."
            );
        } else if (reactants.empty()) {
            addIssue(
                result,
                IssueCategory::constitutiveSource,
                Severity::note,
                reactionIndex,
                "Reaction is a constitutive source and can create material without a represented precursor."
            );
        } else if (products.empty()) {
            addIssue(
                result,
                IssueCategory::terminalSink,
                Severity::note,
                reactionIndex,
                "Reaction is a terminal sink and removes represented material."
            );
        }

        const std::string vectorKey = reactionVectorKey(net);
        if (const auto found = reactionVectors.find(vectorKey); found != reactionVectors.end()) {
            addIssue(
                result,
                IssueCategory::duplicateReactionVector,
                Severity::note,
                reactionIndex,
                "Reaction has the same net stoichiometric vector as reaction " +
                    std::to_string(found->second) + "."
            );
        } else {
            reactionVectors.emplace(vectorKey, reactionIndex);
        }

        for (const auto& [species, coefficient] : net) {
            if ((ir.species[species].flags & speciesConserved) != 0 && coefficient != 0) {
                addIssue(
                    result,
                    IssueCategory::declaredConservationViolation,
                    Severity::error,
                    species,
                    "Species is marked individually conserved but reaction " +
                        std::to_string(reactionIndex) + " changes it."
                );
            }
        }

        if (options.computeComplexDeficiency) {
            const std::uint32_t source = internComplex(canonicalComplex(reactants));
            const std::uint32_t destination = internComplex(canonicalComplex(products));
            result.complexes[source].outgoingComplexes.push_back(destination);
            result.complexes[destination].incomingComplexes.push_back(source);
            complexEdges.emplace_back(source, destination);
        }
    }

    for (std::uint32_t species = 0; species < speciesCount; ++species) {
        if (!consumed[species] && !produced[species]) {
            addIssue(
                result,
                IssueCategory::isolatedSpecies,
                Severity::warning,
                species,
                "Species is not connected to any reaction."
            );
        } else if (!consumed[species]) {
            addIssue(
                result,
                IssueCategory::sourceOnlySpecies,
                Severity::note,
                species,
                "Species is produced but is never consumed by the represented network."
            );
        } else if (!produced[species]) {
            addIssue(
                result,
                IssueCategory::sinkOnlySpecies,
                Severity::note,
                species,
                "Species is consumed but is never produced by the represented network."
            );
        }
    }

    std::map<std::size_t, Component> componentMap;
    for (std::uint32_t species = 0; species < speciesCount; ++species) {
        componentMap[networkComponents.find(species)].speciesIndices.push_back(species);
    }
    for (std::uint32_t reaction = 0; reaction < reactionCount; ++reaction) {
        componentMap[networkComponents.find(speciesCount + reaction)].reactionIndices.push_back(reaction);
    }
    for (auto& [root, component] : componentMap) {
        (void)root;
        result.components.push_back(std::move(component));
    }
    std::sort(result.components.begin(), result.components.end(), [](const Component& left, const Component& right) {
        const auto leftSpecies = left.speciesIndices.empty() ? std::numeric_limits<std::uint32_t>::max() : left.speciesIndices.front();
        const auto rightSpecies = right.speciesIndices.empty() ? std::numeric_limits<std::uint32_t>::max() : right.speciesIndices.front();
        if (leftSpecies != rightSpecies) return leftSpecies < rightSpecies;
        const auto leftReaction = left.reactionIndices.empty() ? std::numeric_limits<std::uint32_t>::max() : left.reactionIndices.front();
        const auto rightReaction = right.reactionIndices.empty() ? std::numeric_limits<std::uint32_t>::max() : right.reactionIndices.front();
        return leftReaction < rightReaction;
    });

    if (options.computeComplexDeficiency && !result.complexes.empty()) {
        result.complexCount = static_cast<std::uint32_t>(result.complexes.size());
        UnionFind linkage(result.complexes.size());
        for (const auto& [source, destination] : complexEdges) {
            linkage.unite(source, destination);
        }
        std::set<std::size_t> linkageRoots;
        for (std::size_t complex = 0; complex < result.complexes.size(); ++complex) {
            linkageRoots.insert(linkage.find(complex));
            auto& record = result.complexes[complex];
            std::sort(record.outgoingComplexes.begin(), record.outgoingComplexes.end());
            record.outgoingComplexes.erase(
                std::unique(record.outgoingComplexes.begin(), record.outgoingComplexes.end()),
                record.outgoingComplexes.end()
            );
            std::sort(record.incomingComplexes.begin(), record.incomingComplexes.end());
            record.incomingComplexes.erase(
                std::unique(record.incomingComplexes.begin(), record.incomingComplexes.end()),
                record.incomingComplexes.end()
            );
        }
        result.linkageClassCount = static_cast<std::uint32_t>(linkageRoots.size());

        std::vector<std::int32_t> indices(result.complexes.size(), -1);
        std::vector<std::int32_t> lowLink(result.complexes.size(), -1);
        std::vector<bool> onStack(result.complexes.size(), false);
        std::vector<std::uint32_t> stack;
        std::vector<std::uint32_t> componentByNode(result.complexes.size(), 0);
        std::int32_t nextIndex = 0;
        std::uint32_t strongCount = 0;
        for (std::uint32_t node = 0; node < result.complexes.size(); ++node) {
            if (indices[node] < 0) {
                tarjanVisit(
                    node,
                    result.complexes,
                    indices,
                    lowLink,
                    onStack,
                    stack,
                    nextIndex,
                    componentByNode,
                    strongCount
                );
            }
        }
        result.strongLinkageClassCount = strongCount;
        result.weaklyReversible = std::all_of(
            complexEdges.begin(),
            complexEdges.end(),
            [&](const auto& edge) {
                return componentByNode[edge.first] == componentByNode[edge.second];
            }
        );
        if (!result.weaklyReversible && !complexEdges.empty()) {
            addIssue(
                result,
                IssueCategory::nonWeaklyReversible,
                Severity::note,
                0,
                "Complex graph is not weakly reversible."
            );
        }
    }

    const std::uint64_t denseElements = static_cast<std::uint64_t>(speciesCount) *
                                        static_cast<std::uint64_t>(reactionCount);
    const bool denseAllowed = speciesCount <= options.maximumDenseSpecies &&
                              reactionCount <= options.maximumDenseReactions &&
                              denseElements <= options.maximumDenseElements;
    if (denseAllowed) {
        std::vector<std::vector<double>> stoichiometry(
            speciesCount,
            std::vector<double>(reactionCount, 0.0)
        );
        for (std::uint32_t reactionIndex = 0; reactionIndex < reactionCount; ++reactionIndex) {
            const ReactionRecord& reaction = ir.reactions[reactionIndex];
            const std::uint64_t reactantEnd = static_cast<std::uint64_t>(reaction.reactantOffset) + reaction.reactantCount;
            const std::uint64_t productEnd = static_cast<std::uint64_t>(reaction.productOffset) + reaction.productCount;
            if (reactantEnd > ir.stoichiometry.size() || productEnd > ir.stoichiometry.size()) continue;
            for (std::uint32_t termIndex = 0; termIndex < reaction.reactantCount; ++termIndex) {
                const auto& term = ir.stoichiometry[reaction.reactantOffset + termIndex];
                if (term.speciesIndex < speciesCount) {
                    stoichiometry[term.speciesIndex][reactionIndex] -= term.coefficient;
                }
            }
            for (std::uint32_t termIndex = 0; termIndex < reaction.productCount; ++termIndex) {
                const auto& term = ir.stoichiometry[reaction.productOffset + termIndex];
                if (term.speciesIndex < speciesCount) {
                    stoichiometry[term.speciesIndex][reactionIndex] += term.coefficient;
                }
            }
        }

        std::vector<std::vector<double>> transpose(
            reactionCount,
            std::vector<double>(speciesCount, 0.0)
        );
        for (std::size_t species = 0; species < speciesCount; ++species) {
            for (std::size_t reaction = 0; reaction < reactionCount; ++reaction) {
                transpose[reaction][species] = stoichiometry[species][reaction];
            }
        }
        const RREFResult decomposition = rref(std::move(transpose), options.rankTolerance);
        result.stoichiometricRank = decomposition.rank;
        result.denseAnalysisComplete = true;
        if (options.inferConservationLaws) {
            result.conservationLaws = inferConservation(
                decomposition,
                stoichiometry,
                options.conservationTolerance,
                options.maximumConservationLaws
            );
        }
        if (decomposition.rank > 0 &&
            decomposition.minimumAcceptedPivot <=
                options.rankTolerance * std::max(1.0, decomposition.maximumInputMagnitude) * 100.0) {
            addIssue(
                result,
                IssueCategory::numericalRankUncertain,
                Severity::warning,
                0,
                "Stoichiometric rank contains a pivot close to the configured tolerance."
            );
        }
    } else {
        addIssue(
            result,
            IssueCategory::denseAnalysisSkipped,
            Severity::note,
            0,
            "Dense rank and conservation analysis was skipped because the bounded matrix budget was exceeded."
        );
    }

    if (options.computeComplexDeficiency && result.denseAnalysisComplete) {
        result.deficiency = static_cast<std::int64_t>(result.complexCount) -
                            static_cast<std::int64_t>(result.linkageClassCount) -
                            static_cast<std::int64_t>(result.stoichiometricRank);
        if (result.deficiency < 0) {
            addIssue(
                result,
                IssueCategory::numericalRankUncertain,
                Severity::warning,
                0,
                "Computed deficiency is negative; topology or numerical rank requires review."
            );
        } else if (result.deficiency > 0) {
            addIssue(
                result,
                IssueCategory::positiveDeficiency,
                Severity::note,
                0,
                "Reaction network has positive structural deficiency " +
                    std::to_string(result.deficiency) + "."
            );
        }
    }

    for (const auto& issue : result.issues) {
        if (issue.severity == Severity::error || issue.severity == Severity::fatal) {
            result.diagnostics.error(
                "NVN100",
                issue.message,
                std::string(issueCategoryName(issue.category)) + ":" +
                    std::to_string(issue.subjectIndex)
            );
        }
    }
    return result;
}

std::string Result::toJson() const {
    std::ostringstream stream;
    stream << "{\"speciesCount\":" << speciesCount
           << ",\"reactionCount\":" << reactionCount
           << ",\"complexCount\":" << complexCount
           << ",\"linkageClassCount\":" << linkageClassCount
           << ",\"strongLinkageClassCount\":" << strongLinkageClassCount
           << ",\"stoichiometricRank\":" << stoichiometricRank
           << ",\"deficiency\":" << deficiency
           << ",\"weaklyReversible\":" << (weaklyReversible ? "true" : "false")
           << ",\"denseAnalysisComplete\":" << (denseAnalysisComplete ? "true" : "false")
           << ",\"components\":[";
    for (std::size_t index = 0; index < components.size(); ++index) {
        if (index != 0) stream << ',';
        stream << "{\"species\":[";
        for (std::size_t item = 0; item < components[index].speciesIndices.size(); ++item) {
            if (item != 0) stream << ',';
            stream << components[index].speciesIndices[item];
        }
        stream << "],\"reactions\":[";
        for (std::size_t item = 0; item < components[index].reactionIndices.size(); ++item) {
            if (item != 0) stream << ',';
            stream << components[index].reactionIndices[item];
        }
        stream << "]}";
    }
    stream << "],\"conservationLaws\":[";
    for (std::size_t lawIndex = 0; lawIndex < conservationLaws.size(); ++lawIndex) {
        if (lawIndex != 0) stream << ',';
        const auto& law = conservationLaws[lawIndex];
        stream << "{\"maximumResidual\":" << law.maximumResidual << ",\"terms\":[";
        for (std::size_t termIndex = 0; termIndex < law.terms.size(); ++termIndex) {
            if (termIndex != 0) stream << ',';
            stream << "{\"speciesIndex\":" << law.terms[termIndex].speciesIndex
                   << ",\"coefficient\":" << law.terms[termIndex].coefficient << '}';
        }
        stream << "]}";
    }
    stream << "],\"issues\":[";
    for (std::size_t index = 0; index < issues.size(); ++index) {
        if (index != 0) stream << ',';
        stream << "{\"category\":";
        writeJsonString(stream, issueCategoryName(issues[index].category));
        stream << ",\"severity\":";
        writeJsonString(stream, severityName(issues[index].severity));
        stream << ",\"subjectIndex\":" << issues[index].subjectIndex
               << ",\"message\":";
        writeJsonString(stream, issues[index].message);
        stream << '}';
    }
    stream << "],\"diagnostics\":" << diagnostics.toJson() << '}';
    return stream.str();
}

std::string_view issueCategoryName(IssueCategory category) noexcept {
    switch (category) {
        case IssueCategory::isolatedSpecies: return "isolated-species";
        case IssueCategory::sourceOnlySpecies: return "source-only-species";
        case IssueCategory::sinkOnlySpecies: return "sink-only-species";
        case IssueCategory::constitutiveSource: return "constitutive-source";
        case IssueCategory::terminalSink: return "terminal-sink";
        case IssueCategory::declaredConservationViolation: return "declared-conservation-violation";
        case IssueCategory::duplicateReactionVector: return "duplicate-reaction-vector";
        case IssueCategory::emptyReaction: return "empty-reaction";
        case IssueCategory::nonWeaklyReversible: return "non-weakly-reversible";
        case IssueCategory::positiveDeficiency: return "positive-deficiency";
        case IssueCategory::denseAnalysisSkipped: return "dense-analysis-skipped";
        case IssueCategory::numericalRankUncertain: return "numerical-rank-uncertain";
        case IssueCategory::invalidTopology: return "invalid-topology";
    }
    return "unknown";
}

} // namespace nvivo::network
