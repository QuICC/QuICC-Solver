/**
 * @file Tools.hpp
 * @brief Implementation of some tools for the profiling timer 
 */

#ifndef QUICC_DEBUG_PROFILER_TOOLS_HPP
#define QUICC_DEBUG_PROFILER_TOOLS_HPP

// Configuration includes
//

// System includes
//
#include <iostream>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Debug {

namespace Profiler {

   /**
    * @brief Implementation of some tools for the profiling timer
    */
   class Tools
   {
      public:
         /**
          * @brief Print profiling output
          */
         static void printInfo(std::ostream* pStream = nullptr);
         
      protected:
         /**
          * @brief Write profiling output
          */
         static void writeTimings(std::ostream& stream);

      private:
         /**
          * @brief Constructor
          */
         Tools();

         /**
          * @brief Destructor
          */
         ~Tools();
   };

}
}
}

#endif // QUICC_DEBUG_PROFILER_TOOLS_HPP
