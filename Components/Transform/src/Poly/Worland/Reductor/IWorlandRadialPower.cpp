/**
 * @file IWorlandRadialPower.cpp
 * @brief Source of the interface to a Worland radial grid power operator (e.g. energy)
 */

// System includes
//

// Project includes
//
#include "QuICC/Transform/Poly/Worland/Reductor/IWorlandRadialPower.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Reduce.hpp"

namespace QuICC {

namespace Transform {

namespace Poly {

namespace Worland {

namespace Reductor {

   IWorlandRadialPower::IWorlandRadialPower()
   {
   }

   void IWorlandRadialPower::initOperators(const Internal::Array& igrid, const Internal::Array& iweights) const
   {
      #if defined QUICC_WORLAND_REDUIMPL_MATRIX
         // Reserve storage for the operators
         this->mOps.reserve(this->mspSetup->slowSize());

         // Loop over harmonic degrees
         for(int i = 0; i < this->mspSetup->slowSize(); i++)
         {
            // Build operator
            this->mOps.push_back(Matrix(this->mspSetup->fastSize(i), igrid.size()));
            Matrix op;
            this->makeOperator(op, igrid, iweights, i);
            this->mOps.back() = op.transpose();
         }

      #elif defined QUICC_WORLAND_REDUIMPL_OTF

         // Store grid and weights
         this->mGrid = igrid;
         this->mWeights = iweights;

      #endif //defined QUICC_WORLAND_REDUIMPL_MATRIX
   }

   void IWorlandRadialPower::applyOperators(Matrix& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<3> fix(this->mProfileTag);

      // assert right sizes for input  matrix
      assert(in.cols() == this->mspSetup->blockSize());
      // assert right sizes for output matrix
      assert(rOut.rows() == this->outRows());
      assert(rOut.cols() == this->outCols());

      int start = 0;
      int outRows = this->mspSetup->fwdSize();
      for(int i = 0; i < this->mspSetup->slowSize(); i++)
      {
         int cols = this->mspSetup->mult(i);
         int inRows = this->mspSetup->fastSize(i);

         this->applyOperator(rOut.block(0,start, outRows, cols), i, in.block(0,start, inRows, cols));

         start += cols;
      }
   }

   void IWorlandRadialPower::applyOperators(MatrixZ& rOut, const MatrixZ& in) const
   {
      throw std::logic_error("Unused interface");
   }

   void IWorlandRadialPower::defaultApplyOperator(Eigen::Ref<Matrix> rOut, const int i, const Eigen::Ref<const MatrixZ>& in) const
   {
      rOut = (this->mOps.at(i).transpose()*in).array().abs2();
   }

   int IWorlandRadialPower::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int IWorlandRadialPower::outCols() const
   {
      return this->mspSetup->blockSize();
   }

} // Reductor
} // Worland
} // Poly
} // Transform
} // QuICC
