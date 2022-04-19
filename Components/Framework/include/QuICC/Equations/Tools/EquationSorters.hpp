/**
 * @file EquationTools.hpp
 * @brief Implementation of equation sorting functors
 */

#ifndef QUICC_EQUATIONS_SORTERS_EQUATIONSORTERS_HPP
#define QUICC_EQUATIONS_SORTERS_EQUATIONSORTERS_HPP

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

namespace Sorters {

   /**
    * @brief Sorting functor by equation type
    */
   class EquationType
   {
      public:
         /**
          * @brief Sort scalar equations by equation type
          *
          * @param eqA Left scalar equation
          * @param eqB Right scalar equation
          */
         bool operator()(std::shared_ptr<EquationData> eqA, std::shared_ptr<EquationData> eqB);

      private:
         /**
          * @brief Convert enum equation type to integer for scalar equation
          *
          * @param eqA Scalar equation
          */
         int computeEquationType(std::shared_ptr<EquationData> eqA);
   };

}
}
}

#endif // QUICC_EQUATIONS_SORTERS_EQUATIONSORTERS_HPP
