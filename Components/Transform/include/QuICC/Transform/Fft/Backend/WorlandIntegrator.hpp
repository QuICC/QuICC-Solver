/**
 * @file WorlandIntegrator.hpp
 * @brief Interface for a generic API for Worland integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_WORLANDINTEGRATOR_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_WORLANDINTEGRATOR_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Transform/Fft/Worland/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   /**
    * @brief Interface for a generic API for Worland integrator
    */
   class WorlandIntegrator
   {
      public:
         /// Typedef for the configuration class
         typedef Worland::Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef Worland::SharedSetup SharedSetupType;

         /**
          * @brief Constructor
          */
         WorlandIntegrator();

         /**
          * @brief Destructor
          */
         virtual ~WorlandIntegrator();

         /**
          * @brief Initialise the FFT transforms
          */
         void init(const SetupType& setup, const int lshift, const bool lshiftOnlyParity = false, const bool alwaysZeroNegative = false) const;

         /**
          * @brief Set zero filter
          */
         void setZFilter(const std::set<int>& filter) const;

         /**
          * @brief Change internal spectral resolution
          */
         void setWSize(const unsigned int shiftMaxL) const;

         /**
          * @brief Initialize addition temporary storage
          */
         void addStorage(const int inExtras, const int outExtras) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         void io(MHDFloat* out, const MHDFloat* in) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         void io(const bool isEven) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         void input(const Matrix& in, const bool isEven) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         void input(const MatrixZ& out, const bool isEven, const bool useReal) const;

         /**
          * @brief Set output
          */
         void output(Matrix& rOut, const bool isEven) const;

         /**
          * @brief Set output
          */
         void output(MatrixZ& rOut, const bool isEven, const bool useReal) const;

         /**
          * @brief Apply FFT
          */
         void applyFft() const;

         /**
          * @brief Convert FFT coefficients to Worland expansion
          */
         void forwardWorland(const bool isEven, const int id = 0) const;

         /**
          * @brief Lower alpha by 1
          */
         void lowerAlpha(const MHDFloat alpha, const bool isEven, const int id = 0, const MHDFloat norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Lower alpha by 1
          */
         void raiseAlpha(const MHDFloat alpha, const bool isEven, const int id = 0, const MHDFloat norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Lower beta by 1
          */
         void lowerBeta(const MHDFloat alpha, const bool isEven, const int id = 0, const MHDFloat norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Raise beta by 1
          */
         void raiseBeta(const MHDFloat alpha, const bool isEven, const int id = 0, const MHDFloat norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Lower beta by 1 and divide by r^2
          */
         void lowerR2Beta(const MHDFloat alpha, const bool isEven, const int id = 0, const MHDFloat norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Raise beta by 1 and multiply by r^2
          */
         void raiseR2Beta(const MHDFloat alpha, const bool isEven, const int id = 0, const MHDFloat norm = 1.0/std::sqrt(2.0)) const;

         /**
          * @brief Apply sparse integration operator I2
          */
         void applyI2(const bool isEven, const int id = 0) const;

         /**
          * @brief Apply sparse integration operator I4
          */
         void applyI4(const bool isEven, const int id = 0) const;

         /**
          * @brief Scaling by constant
          */
         void scaleC(const MHDFloat c, const bool isEven, const unsigned int id = 0) const;

         /**
          * @brief Scale coefficients by a*l + y
          */
         void scaleALPY(const MHDFloat a, const MHDFloat y, const bool isEven, const int lshift = 0, const unsigned int id = 0) const;

         /**
          * @brief Scale coefficients by 2.0(l+n+1)*sqrt((n+1)/(l+n+1))
          */
         void scaleD(const bool isEven, const int lshift = 0, const unsigned int id = 0) const;

         /**
          * @brief Scale coefficients by ??
          */
         void scaleSphLaplA(const bool isEven, const int lshift = 0, const unsigned int id = 0) const;

         /**
          * @brief Scale coefficients by ??
          */
         void scaleSphLaplB(const bool isEven, const int lshift = 0, const unsigned int id = 0) const;

         /**
          * @brief Shift l in temporary storage
          */
         void lshift(const unsigned int id, const int lshift, const bool isEven) const;

         /**
          * @brief Shift n expansion in temporary storage
          */
         void nshift(const unsigned int id, const int nshift, const bool isEven) const;

         /**
          * @brief Copy out temporary storage into new extra temporary
          */
         void copy(const int to, const int from, const int nshift, const bool isEven) const;

         /**
          * @brief Add extra temporary to out temporary storage
          */
         void add(const int to, const int from, const int nshift, const bool isEven) const;

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

#endif // QUICC_TRANSFORM_FFT_BACKEND_WORLANDINTEGRATOR_HPP
