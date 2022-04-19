/**
 * @file Library.cpp
 * @brief Source for the static interface to FFTW
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/Library.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

  // Fastest FFTW plan creation
  #ifdef QUICC_FFTPLAN_FAST
     unsigned int Library::sPlanFlag = FFTW_ESTIMATE;
  // Medium FFTW plan creation
  #elif defined QUICC_FFTPLAN_MEDIUM
     unsigned int Library::sPlanFlag = FFTW_MEASURE;
  // Slow FFTW plan creation
  #elif defined QUICC_FFTPLAN_SLOW
     unsigned int Library::sPlanFlag = FFTW_PATIENT;
  #endif // QUICC_FFTW_ESTIMATE

   int Library::sCounter = 0;

   Library::Library()
   {
   }

   Library::~Library()
   {
   }

   unsigned int Library::planFlag()
   {
      return Library::sPlanFlag;
   }

   void Library::registerUser()
   {
      ++Library::sCounter;
   }

   void Library::unregisterUser()
   {
      --Library::sCounter;
   }

   void Library::init()
   {
      #if defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
      if(Library::sCounter == 0)
      {
         // Initialize FFTW's threads
         int error = fftw_init_threads();

         if(error == 0)
         {
            throw std::logic_error("FFTW's threads initialization failed!");
         }
      }

      // Number of threads
      fftw_plan_with_nthreads(2);
      #endif //defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
   }

   void Library::cleanup()
   {
      // Check if all objects have been destroyed
      if(Library::sCounter == 0)
      {
         #if defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
            fftw_cleanup_threads();
         #else
            fftw_cleanup();
         #endif //defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
      }
   }

}
}
}
}
}
