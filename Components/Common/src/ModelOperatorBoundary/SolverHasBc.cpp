/**
 * @file SolverHasBc.cpp
 * @brief Source of the SolverHasBc ModelOperatorBoundary
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"

// Project includes
//

namespace QuICC {

namespace ModelOperatorBoundary {

   std::string SolverHasBc::sTag()
   {
      return "solver_has_bc";
   }

   std::string SolverHasBc::sFormatted()
   {
      return "SolverHasBC";
   }

   SolverHasBc::SolverHasBc()
      : IRegisterId<SolverHasBc>(SolverHasBc::sTag(), SolverHasBc::sFormatted())
   {
   }

}
}
