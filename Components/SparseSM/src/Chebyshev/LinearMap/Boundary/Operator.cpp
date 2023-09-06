/** 
 * @file Operator.cpp
 * @brief Source of the implementation of boundary operator
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Boundary/Operator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

namespace Boundary {

   Operator::Operator(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const bool atTop)
      : ILinearMapOperator(rows, cols, lower, upper), mAtTop(atTop), mLower(lower), mUpper(upper)
   {
   }

   void Operator::buildTriplets(TripletList_t& list) const
   {
      assert(this->mBcs.size() <= static_cast<std::size_t>(this->rows()));

      int i = 0;
      if(!this->mAtTop)
      {
         i = this->rows()-this->mBcs.size();
      }

      for(auto row: this->mBcs)
      {
         assert(row.size() >= this->cols());

         for(int j = 0; j < this->cols(); j++)
         {
            Triplet_t t(i, j, row(j));
            list.emplace_back(t);
         }
         i++;
      }
   }
 
} // Boundary
} // LinearMap
} // Chebyshev
} // Polynomial
} // QuICC
