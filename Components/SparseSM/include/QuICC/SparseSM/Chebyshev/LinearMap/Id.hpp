/**
 * @file Id.hpp
 * @brief Implementation of the Chebyshev linear map generalized identity sparse operator. Part of the top/bottom can be made zero and it doesn't need to be square.
 */

#ifndef QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_ID_HPP
#define QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_ID_HPP

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

   /**
    * @brief Implementation of the Chebyshev linear map generalized identity sparse operator. Part of the top/bottom can be made zero and it doesn't need to be square.
    */
   class Id: public ILinearMapOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param rows    Number of row
          * @param cols    Number of cols
          * @param lower   Lower bound
          * @param upper   Upper bound
          * @param q       Truncation q (only consider rows - q equations)
          * @param s Shift of main diagonal
          */
         Id(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const int q = 0, const int s = 0);

         /**
          * @brief Destructor
          */
         virtual ~Id() = default;

      protected:

      private:
         /**
          * @brief diagonal
          *
          * @param n mode indexes
          */
         ACoeff_t d0(const ACoeff_t& n) const;

         /**
          * @brief Build triplet representation of matrix
          *
          * @param list List of triplets (row, col, value)
          */
         virtual void buildTriplets(TripletList_t& list) const override;

         /**
          * @brief Truncation Q
          */
         int mQ;

         /**
          * @brief Shift of main diagonal
          */
         int mShift;
   };

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC

#endif // QUICC_SPARSESM_CHEBYSHEV_LINEARMAP_ID_HPP
