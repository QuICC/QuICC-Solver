/**
 * @file EquationOptions.hpp
 * @brief Base class for holding special options for the equation
 */

#ifndef QUICC_EQUATIONS_EQUATIONOPTIONS_HPP
#define QUICC_EQUATIONS_EQUATIONOPTIONS_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Lowest building block for the implementation of an equation
    */
   class EquationOptions
   {
      public:
         /**
          * @brief Default constructor
          */
         EquationOptions();

         /**
          * @brief Simple constructor
          *
          * @param it Iteration index
          */
         explicit EquationOptions(const int it);

         /**
          * @brief Simple empty destructor
          */
         virtual ~EquationOptions() = default;

         /**
          * @brief Sub-iteration
          */
         int it() const;

      protected:

         /**
          * @brief Sub-iteration at which equation is active
          */
        const int mIt;
   };
} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_EQUATIONOPTIONS_HPP
