/**
 * @file Setup.hpp
 * @brief Implementation of the Mixed FFT setup class 
 */

#ifndef QUICC_TRANSFORM_FFT_FOURIER_MIXED_SETUP_HPP
#define QUICC_TRANSFORM_FFT_FOURIER_MIXED_SETUP_HPP

// Configuration includes
//
#include <memory>

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

   /**
    * @brief Implementation of the Mixed FFT setup class
    */ 
   class Setup: public ::QuICC::Transform::Fft::Setup
   {
      public:
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
         
      protected:

      private:
         /**
          * @brief Set backward size
          */
         virtual void setBackwardSize();
   };

   /// Typedef for an smart reference counting pointer for a Setup
   typedef std::shared_ptr<Setup>   SharedSetup;

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_FOURIER_MIXED_SETUP_HPP
