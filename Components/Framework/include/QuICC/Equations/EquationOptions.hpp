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
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base class for holding special option for equation
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
          * @brief Simple constructor
          *
          * @param it         Iteration index
          * @param nlIsLhs    Nonlinear term is LHS?
          * @param traHasQi   Transform includes QI?
          */
         EquationOptions(const int it, const bool nlIsLhs, const bool traHasQi);

         /**
          * @brief Simple constructor
          *
          * @param it         Iteration index
          * @param nlIsLhs    Nonlinear term is LHS?
          * @param traHasQi   Transform includes QI?
          * @param isBase     Is base equation? (ie not Jacobian)
          */
         EquationOptions(const int it, const bool nlIsLhs, const bool traHasQi, const bool isBase);

         /**
          * @brief Simple empty destructor
          */
         virtual ~EquationOptions() = default;

         /**
          * @brief Sub-iteration
          */
         int it() const;

         /**
          * @brief Nonlinear term is on LHS?
          */
        const bool nonlinearIsLhs;

         /**
          * @brief Transform includes QI?
          */
        const bool transformHasQi;

         /**
          * @brief Is base equation? (ie not Jacobian)
          */
        const bool isBase;

      protected:
         /**
          * @brief Sub-iteration at which equation is active
          */
        const int mIt;
   };
} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_EQUATIONOPTIONS_HPP
