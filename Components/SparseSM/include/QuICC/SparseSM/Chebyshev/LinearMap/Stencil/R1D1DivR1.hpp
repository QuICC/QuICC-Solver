/**
 * @file R1D1DivR1.hpp
 * @brief Implementation of the boundary value of r D 1/r stencil
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_STENCIL_R1D1DIVR1_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_STENCIL_R1D1DIVR1_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SparseSM/Chebyshev/ILinearMapOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Stencil {

   /**
    * @brief Implementation of the boundary value of r D 1/r stencil
    */
   class R1D1DivR1: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of rows
          * @param cols    Number of columns
          * @param lower   Lower bound
          * @param upper   Upper bound
          */
         R1D1DivR1(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper);

         /**
          * @brief Destructor
          */
         virtual ~R1D1DivR1() = default;

      protected:

      private:
         /**
          * @brief 2nd subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_2(const ACoeff_t& n) const;

         /**
          * @brief 1st subdiagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d_1(const ACoeff_t& n) const;

         /**
          * @brief diagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d0(const ACoeff_t& n) const;

         /**
          * @brief Build triplet representation of matrix
          *
          * @param[out] list containing triplets
          */
         virtual void buildTriplets(TripletList_t& list) const;
   };

} // Stencil
} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_STENCIL_R1D1DIVR1_HPP
