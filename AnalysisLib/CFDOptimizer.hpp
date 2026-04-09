#pragma once

#include "libCo.hpp"
#include "TH2F.h"
#include "../CaenLib/RootHit.hpp"
#include "TraceAnalysis.hpp"

namespace Caen1725
{
  struct CFDParameters
  {
  public:
    std::vector<double> shifts;
    std::vector<double> fractions;
  };

  class CFDOptimizer
  {
  public:
    CFDOptimizer() noexcept = default;
    
    
  };
}