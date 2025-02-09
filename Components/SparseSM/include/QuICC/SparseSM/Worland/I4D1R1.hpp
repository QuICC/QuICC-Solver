/** 
 * @file I4D1R1.hpp
 * @brief Implementation of the full sphere Worland I4D1R1 sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_I4D1R1_HPP
#define QUICC_SPARSESM_WORLAND_I4D1R1_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/IWorlandOperator.hpp"
#include "QuICC/SparseSM/Worland/I4D1R1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland I4D1R1 sparse operator
    */ 
   class I4D1R1: public IWorlandOperator
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
          * @param q       Truncation q (only consider rows - q equations)
          */
         I4D1R1(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q = 0);

         /**
          * @brief Destructor
          */
         virtual ~I4D1R1() = default;
         
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
         std::shared_ptr<I4D1R1Diags> mpImpl;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_I4D1R1_HPP
