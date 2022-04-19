/**
 * @file Setup.hpp
 * @brief Implementation of the Chebyshev FFT setup class 
 */

#ifndef QUICC_TRANSFORM_FFT_CHEBYSEV_SETUP_HPP
#define QUICC_TRANSFORM_FFT_CHEBYSEV_SETUP_HPP

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

namespace Chebyshev {

   /**
    * @brief Implementation of the Chebyshev FFT setup class
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
          * @brief Set the bounds of the Chebyshev domain y = [lower, upper] wrt x = [-1, 1] (natural Chebyshev grid)
          *
          * @param lower      Lower bound of linear map y = ax + b, x = [-1, 1] (natural Chebyshev grid)
          * @param upper      Upper bound of linear map y = ax + b, x = [-1, 1] (natural Chebyshev grid)
          */
         void setBounds(const MHDFloat lower, const MHDFloat upper);

         /**
          * @brief Lower bound for y = ax + b mapping
          */
         MHDFloat lower() const;

         /**
          * @brief Upper bound for y = ax + b mapping
          */
         MHDFloat upper() const;

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

         /**
          * @brief Lower bound of linear map
          */
         MHDFloat mLower;

         /**
          * @brief Upper bound of linear map
          */
         MHDFloat mUpper;
   };

   /// Typedef for an smart reference counting pointer for a Setup
   typedef std::shared_ptr<Setup>   SharedSetup;

}
}
}
}

#endif // QUICC_TRANSFORM_FFT_CHEBYSHEV_SETUP_HPP
