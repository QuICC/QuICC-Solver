/**
 * @file EquationConditions.hpp
 * @brief Implementation of equation condition functors
 */

#ifndef QUICC_EQUATIONS_CONDITIONS_EQUATIONCONDITIONS_HPP
#define QUICC_EQUATIONS_CONDITIONS_EQUATIONCONDITIONS_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Conditions {

   /**
    * @brief Condition functor for prognostic equation
    */
   struct IsPrognostic
   {
      /**
       * @brief Check scalar equation
       */
      bool operator()(std::shared_ptr<EquationData> eqA);
   };

   /**
    * @brief Condition functor for diagnostic equation
    */
   struct IsDiagnostic
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(std::shared_ptr<EquationData> eqA);
   };

   /**
    * @brief Condition functor for trivial equation
    */
   struct IsTrivial
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(std::shared_ptr<EquationData> eqA);
   };

   /**
    * @brief Condition functor for wrapper
    */
   struct IsWrapper
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(std::shared_ptr<EquationData> eqA);
   };

}
}
}

#endif // QUICC_EQUATIONS_CONDITIONS_EQUATIONCONDITIONS_HPP
