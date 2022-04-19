/**
 * @file CouplingFeature.cpp
 * @brief Source of coupling feature
 */

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/CouplingFeature.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

   std::map<CouplingFeature,bool> defaultCouplingFeature()
   {
      std::map<CouplingFeature,bool> features;
      features.insert(std::make_pair(CouplingFeature::Nonlinear, false));
      features.insert(std::make_pair(CouplingFeature::Source, false));
      features.insert(std::make_pair(CouplingFeature::BoundaryValue, false));
      features.insert(std::make_pair(CouplingFeature::AllowExplicit, true));

      return features;
   }
}
}
