/**
 * @file Setup.hpp
 * @brief Implementation of the FFT setup class 
 */

#ifndef QUICC_TRANSFORM_FFT_SETUP_HPP
#define QUICC_TRANSFORM_FFT_SETUP_HPP

// Configuration includes
//
#include <memory>

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/BasicTypes.hpp"
#include "QuICC/Transform/TransformSetup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

   /**
    * @brief Implementation of the FFT setup class
    */ 
   class Setup: public TransformSetup
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
          * @brief Set the box size scaling factor
          */
         void setBoxScale(const MHDFloat boxScale);

         /**
          * @brief Get the size of the transform
          */
         int bwdSize() const;

         /**
          * @brief Get the size of the padding
          */
         int padSize() const;

         /**
          * @brief Get box size scaling factor
          */
         MHDFloat boxScale() const;
         
      protected:
         /**
          * @brief Size of the backward transform
          */
         int mBwdSize;

      private:
         /**
          * @brief Set backward size
          */
         virtual void setBackwardSize() = 0;

         /**
          * @brief Storage for the box scale factor for derivatives
          */
         MHDFloat mBoxScale;
   };

   /// Typedef for an smart reference counting pointer for a Setup
   typedef std::shared_ptr<Setup>   SharedSetup;

}
}
}

#endif // QUICC_TRANSFORM_FFT_SETUP_HPP
