/**
 * @file Library.cpp
 * @brief Source for the static interface to FFTW
 */

// System includes
//
#include <stdexcept>
#include <fftw3.h>

// Project includes
//
#include "Library.hpp"

namespace QuICC {
namespace Fft {
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

   Library& Library::getInstance()
   {
      static Library instance;
      return instance;
   }

   Library::Library()
   {
      #if defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
      // Initialize FFTW's threads
      int error = fftw_init_threads();

      if(error == 0)
      {
         throw std::logic_error("FFTW's threads initialization failed!");
      }

      // Number of threads
      fftw_plan_with_nthreads(2);
      #endif //defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
   }

   Library::~Library()
   {
      #if defined QUICC_THREADS_PTHREADS || defined QUICC_THREADS_OPENMP
         fftw_cleanup_threads();
      #else
         fftw_cleanup();
      #endif //defined QUICC_THREADS_PTHREADS || defined
   }

   unsigned int Library::planFlag()
   {
      return Library::sPlanFlag;
   }


} // namespace Fftw
} // namespace Fft
} // namespace QuICC
