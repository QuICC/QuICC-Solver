/**
 * @file Setup.hpp
 * @brief Implementation of the Complex FFT setup class
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_SETUP_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_SETUP_HPP

// Configuration includes
//
#include <memory>

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

   /**
    * @brief Implementation of the Complex FFT setup class
    */
   class Setup: public ::QuICC::Transform::Fft::Setup
   {
      public:
         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param specSize   Spectral output size (i.e without the padding)
          */
         Setup(const int size, const int specSize, const GridPurpose::Id purpose);

         /**
          * @brief Constructor
          *
          * @param size       Size of the transform
          * @param blockSize  Number of similar transforms
          * @param specSize   Spectral output size (i.e without the padding)
          */
         Setup(const int size, const int blockSize, const int specSize, const GridPurpose::Id purpose);

         /**
          * @brief Empty destructor
          */
         virtual ~Setup();

         /**
          * @brief Lock setup to forbid adding new indexes
          */
         virtual void lock();

         /**
          * @brief Get the block distribution with 3D index
          */
         const MatrixI& idBlocks() const;

      protected:

      private:
         /**
          * @brief Set backward size
          */
         virtual void setBackwardSize();

         /**
          * @brief Storage for the block sizes with 3D index
          */
         MatrixI mIdBlocks;
   };

   /// Typedef for an smart reference counting pointer for a Setup
   typedef std::shared_ptr<Setup>   SharedSetup;

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_COMPLEX_SETUP_HPP
