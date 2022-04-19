/** 
 * @file ComplexIntegrator.hpp
 * @brief Interface for a generic API for complex integrator 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_COMPLEXINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_COMPLEXINTEGRATOR_HPP

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
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   /**
    * @brief Interface for a generic API for complex integrator
    */ 
   class ComplexIntegrator
   {
      public:
         /// Typedef for the configuration class
         typedef Fourier::Complex::Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef Fourier::Complex::SharedSetup SharedSetupType;

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
         void init(const SetupType& setup) const;

         /**
          * @brief Setup the mean blocks
          */
         void initMeanBlocks(const MatrixI& idBlocks) const;

         /**
          * @brief Set input data pointers for FFT (uses internal pointer for output)
          */
         void input(const MHDComplex* in) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(MHDComplex* out, const MHDComplex* in) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(MatrixZ& rOut, const MatrixZ& in) const;

         /**
          * @brief Copy field out of backend
          */
         void output(MatrixZ& rOut) const;

         /**
          * @brief Set output data pointers for FFT (uses internal pointer for input)
          */
         void output(MHDComplex* out) const;

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
         void outputDiff(MatrixZ& rOut, const int order, const MHDFloat scale) const;

         /**
          * @brief Extract the mean
          */
         void extractMean() const;

         /**
          * @brief Set the mean
          */
         void setMean(MatrixZ& rOut, const MHDFloat scale) const;

         /**
          * @brief Apply FFT
          */
         void applyFft() const;

         /**
          * @brief Compute 2D derivative operator
          */
         int computeDiff2D(const std::vector<std::pair<int,int> >& orders, const MHDFloat scale, const MatrixI& idBlocks, const bool isInverse = false) const;

         /**
          * @brief Multiply two 2D derivative operator
          */
         int multDiff2D(const int idA, const int idB) const;

         /**
          * @brief Compute 2D derivative operator
          */
         void applyDiff2D(MatrixZ& rOut, const int id) const;

         /**
          * @brief Destroy 2D derivative operator
          */
         void destroyDiff2D(const int id) const;

         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
         
      protected:

      private:
         /**
          * @brief PIMPL type forward declaration
          */
         struct BackendImpl;

         /**
          * @brief PIMPL
          */
         std::shared_ptr<BackendImpl> mpImpl;
   };

}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_COMPLEXINTEGRATOR_HPP
