/**
 * @file SolverNoTau.cpp
 * @brief Source of the SolverNoTau ModelOperatorBoundary
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperatorBoundary {

   std::string SolverNoTau::sTag()
   {
      return "solver_no_tau";
   }

   std::string SolverNoTau::sFormatted()
   {
      return "SolverNoTau";
   }

   SolverNoTau::SolverNoTau()
      : IRegisterId<SolverNoTau>(SolverNoTau::sTag(), SolverNoTau::sFormatted())
   {
   }

}
}
