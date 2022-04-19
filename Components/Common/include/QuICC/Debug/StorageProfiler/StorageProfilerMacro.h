/** 
 * @file StorageProfilerMacro.h
 * @brief Preprocessor macros to setup storage profiling if requested in configuration
 */

#ifndef QUICC_DEBUG_STORAGEPROFILERMACRO_H
#define QUICC_DEBUG_STORAGEPROFILERMACRO_H

#ifdef QUICC_STORAGEPROFILE
   // include storage profiler
   #include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
   #include "QuICC/Debug/StorageProfiler/StorageProfiler.hpp"
   #include "QuICC/Debug/StorageProfiler/StorageProfilerTools.hpp"

   namespace QuICC {
      /// Typedef for a profiler
      typedef Debug::StorageProfiler  StorageProfilerMacro;
   }

   /// Define storage profiler update macro function
   #define StorageProfilerMacro_update(P,m)  Debug::StorageProfiler::update(P,m)

   /// Define storage profiler printInfo macro function
   #define StorageProfilerMacro_printInfo()  Debug::StorageProfilerTools::printInfo()

#else
   /// Define storage profiler update macro function
   #define StorageProfilerMacro_update(P,m)  

   /// Define storage profiler printInfo macro function
   #define StorageProfilerMacro_printInfo()  

#endif // QUICC_STORAGEPROFILE

#endif // QUICC_DEBUG_STORAGEPROFILERMACRO_H
