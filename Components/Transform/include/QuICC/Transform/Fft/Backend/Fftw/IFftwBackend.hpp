/** 
 * @file IFftwBackend.hpp
 * @brief Interface for a generic FFTW based backend 
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_IFFTWBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_IFFTWBACKEND_HPP

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
#include "QuICC/Transform/Fft/Backend/Fftw/Library.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic FFTW backend
    */ 
   class IFftwBackend
   {
      public:
         /**
          * @brief Constructor
          */
         IFftwBackend();

         /**
          * @brief Destructor
          */
         virtual ~IFftwBackend();

         /**
          * @brief Apply FFT
          */
         virtual void applyFft() const;
         virtual void applyFft(Matrix&, const MatrixZ&) const; // =0?
         virtual void applyFft(MatrixZ&, const Matrix&) const; // =0?
         virtual void applyFft(MatrixZ&, const MatrixZ&) const; // =0?

         /**
          * @brief Get the memory requirements
          */
         virtual MHDFloat requiredStorage() const;
         
      protected:
         /**
          * @brief Plan for the transform
          */
         mutable fftw_plan   mPlan;

      private:
         /**
          * @brief Initialise the FFTW library
          */
         void initLibrary() const;

         /**
          * @brief Cleanup the FFT
          */
         void cleanupFft();
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_IFFTWBACKEND_HPP
