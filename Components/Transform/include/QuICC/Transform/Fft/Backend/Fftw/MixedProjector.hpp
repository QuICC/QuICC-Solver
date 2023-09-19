/**
 * @file MixedProjector.hpp
 * @brief Interface for a generic mixed FFTW based projector
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDPROJECTOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDPROJECTOR_HPP

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
#include "QuICC/Transform/Fft/Backend/Fftw/IMixedBackend.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic mixed FFTW based projector
    */
   class MixedProjector: public IMixedBackend
   {
      public:
         /**
          * @brief Constructor
          */
         MixedProjector();

         /**
          * @brief Destructor
          */
         ~MixedProjector();

         /**
          * @brief Initialise the FFTW transforms
          */
         void init(const SetupType& setup) const override;

         /**
          * @brief Copy and pad FFT temporary input
          */
         void input(MatrixZ& out, const MatrixZ& in) const;

         /**
          * @brief Apply spectral derivative of given order
          * and and pad FFT temporary input
          */
         void inputDiff(MatrixZ& out, const MatrixZ& rData, const int order, const MHDFloat scale) const;

         /**
          * @brief Apply FFT
          */
         void applyFft(Matrix& phys, const MatrixZ& mods) const override;

      protected:

      private:
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

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_MIXEDPROJECTOR_HPP
