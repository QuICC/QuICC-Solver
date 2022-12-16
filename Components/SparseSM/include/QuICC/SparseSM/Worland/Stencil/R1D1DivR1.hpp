/**
 * @file R1D1DivR1.hpp
 * @brief Implementation of the full sphere Worland R1D1DivR1 boundary condition stencil
 */

#ifndef QUICC_SPARSESM_WORLAND_STENCIL_R1D1DIVR1_HPP
#define QUICC_SPARSESM_WORLAND_STENCIL_R1D1DIVR1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/IWorlandOperator.hpp"
#include "QuICC/SparseSM/Worland/Stencil/R1D1DivR1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   /**
    * @brief Implementation of the full sphere Worland R1D1DivR1 boundary condition stencil
    */
   class R1D1DivR1: public IWorlandOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param alpha   Jacobi alpha
          * @param dBeta   Jacobi beta = l + dBeta
          * @param l       Harmonic degree l
          */
         R1D1DivR1(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l);

         /**
          * @brief Destructor
          */
         virtual ~R1D1DivR1() = default;

      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          *
          * @param list List of triplets (row, col, value)
          */
         void buildTriplets(TripletList_t& list) const final;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<R1D1DivR1Diags> mpImpl;
   };

} // Stencil
} // Worland
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_WORLAND_STENCIL_R1D1DIVR1_HPP
