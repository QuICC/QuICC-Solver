/** 
 * @file ImExSBDF2.hpp
 * @brief Implementation of an implicit/explicit SBDF scheme of order 2
 */

#ifndef QUICC_TIMESTEP_IMEXSBDF2_HPP
#define QUICC_TIMESTEP_IMEXSBDF2_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Timestep/SparseOldImExTimestepper.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of an implicit/explicit SBDF scheme of order 2
    */
   class ImExSBDF2: public IImExOldScheme
   {
      public:
         /// Typedef for underlying scheme implementation
         template <typename T1,typename T2,template <typename> class T3>
            using ImplementationType = SparseOldImExTimestepper<T1,T2,T3>;

         /**
          * @brief Constructor
          */
         ImExSBDF2();

         /**
          * @brief Destructor
          */
         virtual ~ImExSBDF2();

         /**
          * @brief Get factor for mass matrix on LHS
          *
          * @param step Substep
          */
         MHDFloat lhsT(const int step) const;

         /**
          * @brief Get factor for linear matrix on LHS
          *
          * @param step Substep
          */
         MHDFloat lhsL(const int step) const;

         /**
          * @brief Get factor for mass matrix on RHS at t_(n-i)
          *
          * @param time Time index
          * @param step Substep
          */
         MHDFloat rhsT(const int i, const int step) const;

         /**
          * @brief Get factor for linear matrix on RHS at t_(n-i)
          *
          * @param time Time index
          * @param step Substep
          */
         MHDFloat rhsL(const int i, const int step) const;

         /**
          * @brief Get factor for nonlinear term at t_(n-i)
          *
          * @param i    Time index
          * @param step Substep
          */
         MHDFloat rhsN(const int i, const int step) const;

         /**
          * @brief Field memory needed at current step
          *
          * @param step Substep
          */
         int fieldMemory(const int step) const;

         /**
          * @brief Nonlinear term memory needed at current step
          *
          * @param step Substep
          */
         int nonlinearMemory(const int step) const;

         /**
          * @brief Butcher's tableau c_i factor for explicit scheme
          */
         MHDFloat cEx(const int i) const;

         /**
          * @brief Number of substeps for final step
          */
         int steps() const;

         /**
          * @brief Order of the scheme
          */
         int order() const;

         /**
          * @brief Scheme has embedded lower order scheme?
          */
         bool hasEmbedded() const;

         /**
          * @brief Number of previous field values required
          */
         int fieldMemory() const;

         /**
          * @brief Number of previous nonlinear terms required
          */
         int nonlinearMemory() const;

         /**
          * @brief Name of the scheme
          */
         std::string name() const;

         /**
          * @brief Initialize Butcher's tableau
          */
         void init();
         
      protected:

      private:

   };

   /// Typedef for a shared pointer ImExSBDF2
   typedef std::shared_ptr<ImExSBDF2> SharedImExSBDF2;
}
}

#endif // QUICC_TIMESTEP_IMEXSBDF2_HPP
