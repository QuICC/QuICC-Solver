/** 
 * @file ComplexIntegrator.hpp
 * @brief Interface for a generic complex cuFFT based integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_COMPLEXINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_COMPLEXINTEGRATOR_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/IComplexBackend.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic complex cuFFT based integrator
    */ 
   class ComplexIntegrator: public IComplexBackend
   {
      public:
         /**
          * @brief Constructor
          */
         ComplexIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~ComplexIntegrator();
         
         /**
          * @brief Initialise the FFT transforms
          */
         virtual void init(const SetupType& setup) const override;

         /**
          * @brief Set input and output data pointers for FFT
          */
         virtual void io(MHDComplex* out, const MHDComplex* in) const override;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(MatrixZ& rOut, const MatrixZ& in) const;

         /**
          * @brief Copy field out of backend
          */
         void output(MatrixZ& rOut) const;
         using IComplexBackend::output;

         /**
          * @brief Copy mean of field out of backend
          */
         void outputMean(MatrixZ& rOut) const;

         /**
          * @brief Copy field out of backend zeroing the mean
          */
         void outputZeroMean(MatrixZ& rOut) const;

         /**
          * @brief Scale with fast index dependent function
          */
         void outputDiff(MatrixZ& rOut, const int order, const double scale) const;

         /**
          * @brief Extract the mean
          */
         void extractMean() const;

         /**
          * @brief Set the mean
          */
         void setMean(MatrixZ& rOut, const double scale) const;

         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const override;

         /**
          * @brief Compute 2D derivative operator
          */
         int computeDiff2D(const std::vector<std::pair<int,int> >& orders, const double scale, const MatrixI& idBlocks, const bool isInverse = false) const;

         /**
          * @brief Multiply two 2D derivative operator
          */
         int multDiff2D(const int idA, const int idB) const;

         /**
          * @brief Apply 2D derivative operator
          */
         void applyDiff2D(MatrixZ& rOut, const int id) const;

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
          * @brief FFT scaling factor
          */
         mutable double mFftScaling;

         /**
          * @brief Next Diff2DOp id
          */
         mutable int mNextDiff2DId;

         /**
          * @brief Map for output data
          */
         mutable Eigen::Map<MatrixZ> mOutMap;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_COMPLEXINTEGRATOR_HPP
