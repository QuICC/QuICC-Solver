/** 
 * @file ComplexProjector.hpp
 * @brief Backend for a generic API for complex projector 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_COMPLEXPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_COMPLEXPROJECTOR_HPP

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
    * @brief Backend for a generic API for complex projector
    */ 
   class ComplexProjector
   {
      public:
         /// Typedef for the configuration class
         typedef Fourier::Complex::Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef Fourier::Complex::SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         ComplexProjector();

         /**
          * @brief Destructor
          */
         virtual ~ComplexProjector();
         
         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup) const;

         /**
          * @brief Setup the mean blocks
          */
         void initMeanBlocks(const MatrixI& idBlocks) const;

         /**
          * @brief Copy field into backend
          */
         void input(const MatrixZ& in) const;

         /**
          * @brief Set input data pointers for FFT (uses internal pointer for output)
          */
         void input(const MHDComplex* in) const;

         /**
          * @brief Copy mean of field into backend
          */
         void inputMean(const MatrixZ& in) const;

         /**
          * @brief Set output data pointers for FFT (uses internal pointer for input)
          */
         void output(MHDComplex* out) const;

         /**
          * @brief Set input and output data pointers for FFT
          */
         void io(MHDComplex* out, const MHDComplex* in) const;

         /**
          * @brief Scale with fast index dependent function
          */
         void inputDiff(const MatrixZ& rData, const int order, const MHDFloat scale) const;

         /**
          * @brief Scale with fast and slow index dependent functin
          */
         void inputDiff2D(const MatrixZ& rData, const std::vector<std::pair<int,int> >& orders, const MHDFloat scale, const MatrixI& idBlocks) const;

         /**
          * @brief Apply FFT
          */
         void applyFft() const;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_COMPLEXPROJECTOR_HPP
