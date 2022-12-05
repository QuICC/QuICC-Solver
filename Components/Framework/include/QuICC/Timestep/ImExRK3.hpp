/**
 * @file ImExRK3.hpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order ~3
 */

#ifndef QUICC_TIMESTEP_IMEXRK3_HPP
#define QUICC_TIMESTEP_IMEXRK3_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Timestep/IImExOldScheme.hpp"
#include "QuICC/Timestep/SparseOldImExTimestepper.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order ~3
    */
   class ImExRK3: public IImExOldScheme
   {
      public:
         /// Typedef for underlying scheme implementation
         template <typename T1,typename T2,template <typename> class T3>
            using ImplementationType = SparseOldImExTimestepper<T1,T2,T3>;

         /**
          * @brief Constructor
          */
         ImExRK3();

         /**
          * @brief Destructor
          */
         virtual ~ImExRK3() = default;

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

      protected:
         /**
          * @brief Initialize coefficients
          */
         void init();

         /**
          * @brief Storage for the alpha parameters
          */
         const Eigen::Array<MHDFloat,3,1> mAlpha;

         /**
          * @brief Storage for the beta parameters
          */
         const Eigen::Array<MHDFloat,3,1> mBeta;

         /**
          * @brief Storage for the gamma parameters
          */
         const Eigen::Array<MHDFloat,3,1> mGamma;

         /**
          * @brief Storage for the zeta parameters
          */
         const Eigen::Array<MHDFloat,3,1> mZeta;

         /**
          * @brief Storage for the explicit c factors
          */
         const Eigen::Array<MHDFloat,3,1> mCEx;

      private:

   };

   /// Typedef for a shared pointer ImExRKCB3
   typedef std::shared_ptr<ImExRK3> SharedImExRK3;
}
}

#endif // QUICC_TIMESTEP_IMEXRK3_HPP
