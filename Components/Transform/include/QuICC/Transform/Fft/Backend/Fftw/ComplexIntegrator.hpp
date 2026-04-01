/**
 * @file ComplexIntegrator.hpp
 * @brief Interface for a generic complex FFTW based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_COMPLEXINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_COMPLEXINTEGRATOR_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/IComplexBackend.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic complex FFTW based integrator
    */
   class ComplexIntegrator: public IComplexBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ComplexIntegrator() = default;

         /**
          * @brief Destructor
          */
         ~ComplexIntegrator() = default;

         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup) const override;

         /**
          * @brief Scale field
          */
         void output(Eigen::Ref<MatrixZ> rOut) const;

         /**
          * @brief Copy mean of field out of backend
          *
          * @param rOut dims: N_k1 x (N_k2 x N_r)
          */
         void outputMean(Eigen::Ref<MatrixZ> rOut) const;

         /**
          * @brief Zero the mean
          */
         void zeroMean(Eigen::Ref<MatrixZ> rOut) const;

         /**
          * @brief Scale with fast index dependent function
          */
         void outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale) const;

         /**
          * @brief Extract the mean
          */
         void extractMean(const MatrixZ& rOut) const;

         /**
          * @brief Set the mean
          */
         void setMean(Eigen::Ref<MatrixZ> rOut, const MHDFloat scale) const;

         /**
          * @brief Compute 2D derivative operator
          */
         int computeDiff2D(const std::vector<std::pair<int,int> >& orders, const MHDFloat scale, const MatrixI& idBlocks, const bool isInverse = false) const;

         /**
          * @brief Multiply two 2D derivative operator
          */
         int multDiff2D(const int idA, const int idB) const;

         /**
          * @brief Apply 2D derivative operator
          */
         void applyDiff2D(Eigen::Ref<MatrixZ> rOut, const int id) const;

         /**
          * @brief Destroy 2D derivative operator
          */
         void destroyDiff2D(const int id) const;

      protected:
         /**
          * brief Initialize new 2D derivative operator
          */
         int addDiff2DOp(const MatrixI& idBlocks) const;

      private:
         typedef std::tuple<MatrixZ,MatrixZ,MatrixI> Diff2DOpType;
         typedef std::map<int,Diff2DOpType> StoreDiff2DOpType;

         /**
          * @brief Bwd size
          */
         mutable int mBwdSize;

         /**
          * @brief Block size
          */
         mutable int mBlockSize;

         /**
          * @brief FFT scaling factor
          */
         mutable MHDFloat mFftScaling;

         /**
          * @brief Next Diff2DOp id
          */
         mutable int mNextDiff2DId;

         /**
          * @brief Operators for 2D derivatives
          */
         mutable StoreDiff2DOpType mDiff2DOp;

         /**
          * @brief Temporary storage for extracted mean
          */
         mutable std::vector<ArrayZ> mTmpMean;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_COMPLEXINTEGRATOR_HPP
