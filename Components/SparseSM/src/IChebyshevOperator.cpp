/** 
 * @file IChebyshevOperator.cpp
 * @brief Source of the implementation of generic interface to a Chebyshev based sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//
#include <unsupported/Eigen/SpecialFunctions>

// Class include
//
#include "QuICC/SparseSM/IChebyshevOperator.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace SparseSM {

   IChebyshevOperator::IChebyshevOperator(const int rows, const int cols)
      : ISparseSMOperator(rows, cols)
   {
   }

   IChebyshevOperator::~IChebyshevOperator()
   {
   }

   void IChebyshevOperator::leftOutOfMatrix(TripletList_t& list, const int row, const int col, const Scalar_t value) const
   {
      list.push_back(Triplet_t(row, std::abs(col), value));
   }

}
}
