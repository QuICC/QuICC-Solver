/**
 * @file Library.hpp
 * @brief Static interface to the global features of the FFTW library
 */
#ifndef QUICC_FFT_FFTW_LIBRARY_HPP
#define QUICC_FFT_FFTW_LIBRARY_HPP

// External includes
//

// Project includes
//

namespace QuICC {
namespace Fft {
namespace Fftw {
   /**
    * @brief Static interface to the global features of the FFTW library
    */
   class Library
   {
      public:
         Library(Library const&) = delete;
         void operator=(Library const&) = delete;

         /**
          * @brief Return singleton
          */
         static Library& getInstance();

         /**
          * @brief Get the plan flag
          */
         static unsigned int planFlag();

      private:

         /**
          * @brief FFTW3 flags for the plan setup
          */
         static unsigned int  sPlanFlag;

         /**
          * @brief Empty constructor
          */
         Library();

         /**
          * @brief Empty destructor
          */
         ~Library();
   };

} // namespace Fftw
} // namespace Fft
} // namespace QuICC

#endif // QUICC_FFT_FFTW_LIBRARY_HPP
