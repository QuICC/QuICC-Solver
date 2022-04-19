/**
 * @file Setup.hpp
 * @brief Implementation of the Worland FFT setup class 
 */

#ifndef QUICC_TRANSFORM_FFT_WORLAND_SETUP_HPP
#define QUICC_TRANSFORM_FFT_WORLAND_SETUP_HPP

// Configuration includes
//
#include <memory>

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Transform/Fft/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

   /**
    * @brief Implementation of the Worland FFT setup class
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

#endif // QUICC_TRANSFORM_FFT_WORLAND_SETUP_HPP
