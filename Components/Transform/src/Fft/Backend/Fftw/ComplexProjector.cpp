/**
 * @file ComplexProjector.cpp
 * @brief Source of the interface for a generic FFTW based complex projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/ComplexProjector.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   ComplexProjector::ComplexProjector()
   {
   }

   ComplexProjector::~ComplexProjector()
   {
   }

   void ComplexProjector::init(const SetupType& setup) const
   {
      //Initialize parent
      IComplexBackend::init(setup);

      int fwdSize = setup.fwdSize();
      int bwdSize = setup.bwdSize();
      int blockSize = setup.blockSize();

      // Initialise temporary storage
      this->mTmp.setZero(bwdSize, blockSize);

      this->mPadSize = setup.padSize();

      this->initMeanBlocks(setup.idBlocks());

      ////////////////////////////////////////
      // Create the plan
      const int  *fftSize = &fwdSize;
      MatrixZ   tmpCplxA(fwdSize, blockSize);
      MatrixZ   tmpCplxB(bwdSize, blockSize);

      this->mPlan = fftw_plan_many_dft(1, fftSize, blockSize, reinterpret_cast<fftw_complex* >(tmpCplxB.data()), NULL, 1, bwdSize, reinterpret_cast<fftw_complex* >(tmpCplxA.data()), NULL, 1, fwdSize, FFTW_BACKWARD, Library::planFlag() | FFTW_DESTROY_INPUT);
      if(this->mPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }
   }

   void ComplexProjector::applyFft(MatrixZ& phys, const MatrixZ& mods) const
   {
      fftw_execute_dft(this->mPlan, reinterpret_cast<fftw_complex *>(const_cast<MHDComplex*>(mods.data())), reinterpret_cast<fftw_complex *>(phys.data()));
   }

   void ComplexProjector::input(MatrixZ& tmp, const MatrixZ& in) const
   {
      tmp.topRows(this->mPosN) = in.topRows(this->mPosN);
      tmp.bottomRows(this->mNegN) = in.bottomRows(this->mNegN);

      this->forceConjugate(tmp);
      this->applyPadding(tmp);
   }

   void ComplexProjector::inputMean(MatrixZ& tmp, const MatrixZ& in) const
   {
      // Initialize to zero
      tmp.setZero();

      // Set the mean
      for(auto it = this->mMeanBlocks.cbegin(); it != this->mMeanBlocks.cend(); ++it)
      {
         tmp.block(0, it->first, 1, it->second) = in.block(0, it->first, 1, it->second);
      }
   }

   void ComplexProjector::inputDiff(MatrixZ& tmp, const MatrixZ& in, const int order, const MHDFloat scale) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*Math::cI;
         tmp.topRows(this->mPosN) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mPosN);
         tmp.bottomRows(this->mNegN) = (sgn*(scale*this->negativeK()).array().pow(order).matrix()).asDiagonal()*in.bottomRows(this->mNegN);
      } else
      {
         MHDFloat sgn = std::pow(-1.0,(order/2)%2);
         tmp.topRows(this->mPosN) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mPosN);
         tmp.bottomRows(this->mNegN) = (sgn*(scale*this->negativeK()).array().pow(order).matrix()).asDiagonal()*in.bottomRows(this->mNegN);
      }

      this->forceConjugate(tmp);
      this->applyPadding(tmp);
   }

   void ComplexProjector::inputDiff2D(MatrixZ& tmp, const MatrixZ& in, const std::vector<std::pair<int,int> >& orders, const MHDFloat scale, const MatrixI& idBlocks) const
   {
      assert(idBlocks.rows() > 0);
      assert(this->mPosN + this->mNegN <= in.rows());

      tmp.topRows(this->mPosN).setZero();
      tmp.bottomRows(this->mNegN).setZero();

      bool isComplex;
      MHDFloat sgn;
      Array fp(this->mPosN);
      Array fn(this->mNegN);
      ArrayZ sZp;
      ArrayZ sZn;
      for(auto& o: orders)
      {
         int fO = o.first;
         int sO = o.second;
         if(fO > 0)
         {
            if(fO%2 == 0)
            {
               sgn = std::pow(-1.0,(fO/2)%2);
               isComplex = false;
            } else
            {
               sgn = std::pow(-1.0,((fO-1)/2)%2);
               isComplex = true;
            }
            fp = (sgn*(scale*this->positiveK()).array().pow(fO).matrix());
            fn = (sgn*(scale*this->negativeK()).array().pow(fO).matrix());
         } else
         {
            fp.setOnes(this->mPosN);
            fn.setOnes(this->mNegN);
            isComplex = false;
         }

         if(sO%2 == 0)
         {
            sgn = std::pow(-1.0,(sO/2)%2);
         } else
         {
            sgn = std::pow(-1.0,((sO-1)/2)%2);
         }
         isComplex = isComplex ^ (sO%2 == 1);

         if(isComplex)
         {
            sZp.resize(fp.size());
            sZp.real().setZero();
            sZn.resize(fn.size());
            sZn.real().setZero();
         }

         int start = 0;
         int negRow = this->mTmp.rows() - this->mNegN;
         for(int i = 0; i < idBlocks.rows(); ++i)
         {
            MHDFloat k = sgn*std::pow(scale*idBlocks(i,0),sO);
            if((fO + sO)%2 == 1)
            {
               sZp.imag() = k*fp;
               sZn.imag() = k*fn;
            } else if(fO%2 == 1)
            {
               k = -k;
            }

            // Split positive and negative frequencies
            if(isComplex)
            {
               tmp.block(0, start, this->mPosN, idBlocks(i,1)) += sZp.asDiagonal()*in.block(0, start, this->mPosN, idBlocks(i,1));
               tmp.block(negRow, start, this->mNegN, idBlocks(i,1)) += sZn.asDiagonal()*in.block(this->mPosN, start, this->mNegN, idBlocks(i,1));
            } else
            {
               tmp.block(0, start, this->mPosN, idBlocks(i,1)) += (k*fp.array()).matrix().asDiagonal()*in.block(0, start, this->mPosN, idBlocks(i,1));
               tmp.block(negRow, start, this->mNegN, idBlocks(i,1)) += (k*fn.array()).matrix().asDiagonal()*in.block(this->mPosN, start, this->mNegN, idBlocks(i,1));
            }

            // Increment block counter
            start += idBlocks(i,1);
         }
      }

      this->forceConjugate(tmp);
      this->applyPadding(tmp);
   }

   void ComplexProjector::forceConjugate(MatrixZ& rData) const
   {
      int endN = rData.rows();

      // Loop over mean blocks
      for(auto it = this->mMeanBlocks.cbegin(); it != this->mMeanBlocks.cend(); ++it)
      {
         // Copy complex conjugate into negative frequency part
         for(int i = 1; i < this->mPosN; i++)
         {
            rData.block(endN - i, it->first, 1, it->second) = rData.block(i, it->first, 1, it->second).conjugate();
         }
      }
   }

   void ComplexProjector::applyPadding(MatrixZ& rData) const
   {
      rData.block(this->mPosN, 0, this->mPadSize, rData.cols()).setZero();
   }

}
}
}
}
}
