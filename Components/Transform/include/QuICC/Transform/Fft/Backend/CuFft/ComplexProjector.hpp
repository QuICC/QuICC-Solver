/**
 * @file ComplexProjector.hpp
 * @brief Backend for a generic complex cuFFT based projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_COMPLEXPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_COMPLEXPROJECTOR_HPP

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
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/IComplexBackend.hpp"
#include "QuICC/Transform/Fft/Fourier/Complex/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Backend for a generic complex cuFFT based projector
    */
   class ComplexProjector: public IComplexBackend
   {
      public:
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
         virtual void init(const SetupType& setup) const override;

         /**
          * @brief Copy field into backend
          */
         void input(const MatrixZ& in) const;
         using IComplexBackend::input;

         /**
          * @brief Copy mean of field into backend
          */
         void inputMean(const MatrixZ& in) const;

         /**
          * @brief Scale with fast index dependent function
          */
         void inputDiff(const MatrixZ& rData, const int order, const double scale) const;

         /**
          * @brief Scale with fast and slow index dependent functin
          */
         void inputDiff2D(const MatrixZ& rData, const std::vector<std::pair<int,int> >& orders, const double scale, const MatrixI& idBlocks) const;

         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const override;

      protected:

      private:
         /**
          * @brief Force exact complex conjugate for mean
          */
         void forceConjugate(MatrixZ& rData) const;

         /**
          * @brief Apply padding
          */
         void applyPadding(MatrixZ& rData) const;

         /**
          * @brief Padding size
          */
         mutable int mPadSize;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_COMPLEXPROJECTOR_HPP
