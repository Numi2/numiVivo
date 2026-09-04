#pragma once

#include "NumiVivoCore/ReactionNetworkAnalysis.hpp"

#include <iterator>

namespace nvivo::network {

class ScalableAnalyzer {
public:
    [[nodiscard]] Result analyze(const ProgramIR& ir,
                                 const Options& options = {}) const;
};

} // namespace nvivo::network
