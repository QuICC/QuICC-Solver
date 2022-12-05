/**
 * @file IImExOldScheme.hpp
 * @brief Interface for a generic implicit/explicit scheme based on old setup
 */

#ifndef QUICC_TIMESTEP_IIMEXOLDSCHEME_HPP
#define QUICC_TIMESTEP_IIMEXOLDSCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Timestep/IScheme.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Interface of generic of an implicit/explicit scheme based on old setup
    */
   class IImExOldScheme: public IScheme
   {
      public:
         /**
          * @brief Constructor
          */
         IImExOldScheme();

         /**
          * @brief Destructor
          */
         virtual ~IImExOldScheme() = default;

         /**
          * @brief Get factor for mass matrix on LHS
          *
          * @param step Substep
          */
         virtual MHDFloat lhsT(const int step) const = 0;

         /**
          * @brief Get factor for linear matrix on LHS
          *
          * @param step Substep
          */
         virtual MHDFloat lhsL(const int step) const = 0;

         /**
          * @brief Get factor for mass matrix on RHS at t_(n-i)
          *
          * @param time Time index
          * @param step Substep
          */
         virtual MHDFloat rhsT(const int i, const int step) const = 0;

         /**
          * @brief Get factor for linear matrix on RHS at t_(n-i)
          *
          * @param time Time index
          * @param step Substep
          */
         virtual MHDFloat rhsL(const int i, const int step) const = 0;

         /**
          * @brief Get factor for nonlinear term at t_(n-i)
          *
          * @param i    Time index
          * @param step Substep
          */
         virtual MHDFloat rhsN(const int i, const int step) const = 0;

         /**
          * @brief Field memory needed at current step
          *
          * @param step Substep
          */
         virtual int fieldMemory(const int step) const = 0;

         /**
          * @brief Nonlinear term memory needed at current step
          *
          * @param step Substep
          */
         virtual int nonlinearMemory(const int step) const = 0;

         /**
          * @brief Butcher's tableau c_i factor for explicit scheme
          */
         virtual MHDFloat cEx(const int i) const = 0;

         /**
          * @brief Number of previous field values required
          */
         virtual int fieldMemory() const = 0;

         /**
          * @brief Number of previous nonlinear terms required
          */
         virtual int nonlinearMemory() const = 0;

      protected:

      private:

   };

   /// Typedef for a shared pointer IImExOldScheme
   typedef std::shared_ptr<IImExOldScheme> SharedIImExOldScheme;

}
}

#endif // QUICC_TIMESTEP_IIMEXOLDSCHEME_HPP
