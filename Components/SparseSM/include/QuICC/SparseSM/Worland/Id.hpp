/**
 * @file Id.hpp
 * @brief Implementation of the full sphere Worland (restricted) identity sparse operator
 */

#ifndef QUICC_SPARSESM_WORLAND_ID_HPP
#define QUICC_SPARSESM_WORLAND_ID_HPP

// System includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SparseSM/IWorlandOperator.hpp"
#include "QuICC/SparseSM/Worland/IdDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   /**
    * @brief Implementation of the full sphere Worland (restricted) identity sparse operator
    */
   class Id: public IWorlandOperator
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
          * @param s Shift of main diagonal
          */
         Id(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q = 0, const int s = 0);

         /**
          * @brief Destructor
          */
         virtual ~Id() = default;

      protected:

      private:
         /**
          * @brief Build triplet representation of matrix
          *
          * @param list List of triplets (row, col, value)
          */
         virtual void buildTriplets(TripletList_t& list) const override;

         /**
          * @brief Implementation of the diagonals
          */
         std::shared_ptr<IdDiags> mpImpl;

         /**
          * @brief Shift of main diagonal
          */
         int mShift;
   };

}
}
}

#endif // QUICC_SPARSESM_WORLAND_ID_HPP
