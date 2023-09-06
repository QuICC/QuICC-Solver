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
#include <fftw3.h>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "Fft/Fftw/Library.hpp"

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
          *
          * deprecated
          */
         virtual void applyFft() const;

         /**
          * @brief Apply FFT
          *
          * Real to real
          */
         virtual void applyFft(Matrix&, const Matrix&) const; //

         /**
          * @brief Apply FFT
          *
          * Complex to real
          */
         virtual void applyFft(Matrix&, const MatrixZ&) const; // =0?

         /**
          * @brief Apply FFT
          *
          * Real to complex
          */
         virtual void applyFft(MatrixZ&, const Matrix&) const; // =0?

         /**
          * @brief Apply FFT
          *
          * Complex to complex
          */
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
