/** 
 * @file IWorlandOperator.cpp
 * @brief Source of the implementation of generic interface to a full sphere Worland sparse operator
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
#include "QuICC/SparseSM/IWorlandOperator.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace SparseSM {

   IWorlandOperator::IWorlandOperator(const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta)
      : ISparseSMOperator(rows, cols)
   {
      if(alpha == MHD_MP(-0.5))
      {
         this->mType = CHEBYSHEV;
      } else if(alpha == MHD_MP(0.0))
      {
         if(dBeta == MHD_MP(-0.5))
         {
            this->mType = LEGENDRE;
         } else if(dBeta == MHD_MP(0.0))
         {
            this->mType = CYLENERGY;
         } else if(dBeta == MHD_MP(0.5))
         {
            this->mType = SPHENERGY;
         } else
         {
            throw std::logic_error("Unknown Legendre type Worland for SparseSM");
         }
      } else
      {
         throw std::logic_error("Unknown Worland type for SparseSM");
      }
   }

   IWorlandOperator::~IWorlandOperator()
   {
   }

   IWorlandOperator::WorlandType IWorlandOperator::type() const
   {
      return this->mType;
   }

}
}
