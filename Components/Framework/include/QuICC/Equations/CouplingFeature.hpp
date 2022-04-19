/**
 * @file CouplingFeature.hpp
 * @brief Enum to some features
 */

#ifndef QUICC_EQUATIONS_COUPLINGFEATURE_HPP
#define QUICC_EQUATIONS_COUPLINGFEATURE_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Equations {

   /**
    * @brief Enum to some features
    */
   enum class CouplingFeature {
      /// Equation has nonlinear terms
      Nonlinear,
      /// Quasi-inverse should be applied
      Source,
      /// Has a boundary value
      BoundaryValue,
      /// Allows explicit terms
      AllowExplicit
   };

   std::map<CouplingFeature,bool> defaultCouplingFeature();

}
}

#endif // QUICC_EQUATIONS_COUPLINGFEATURE_HPP
