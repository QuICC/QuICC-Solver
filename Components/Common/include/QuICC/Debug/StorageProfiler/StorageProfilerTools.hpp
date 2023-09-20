/**
 * @file StorageProfilerTools.hpp
 * @brief Implementation of some tools for the storage profiler
 */

#ifndef QUICC_DEBUG_STORAGEPROFILERTOOLS_HPP
#define QUICC_DEBUG_STORAGEPROFILERTOOLS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Debug {

   /**
    * @brief Implementation of some toos for the strorage profiler
    */
   class StorageProfilerTools
   {
      public:
         /**
          * @brief Print profiling output
          */
         static void printInfo();

         /**
          * @brief Define the unit of the memory requirements
          */
         static void setUnit(MHDFloat value, std::string &ext, MHDFloat &factor);

      protected:

      private:
         /**
          * @brief Constructor
          */
         StorageProfilerTools();

         /**
          * @brief Destructor
          */
         ~StorageProfilerTools();
   };

}
}

#endif // QUICC_DEBUG_STORAGEPROFILERTOOLS_HPP
