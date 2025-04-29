/**
 * @file StorageProfilerTools.cpp
 * @brief Source of the implementation of tools for the storage profiler
 */

// Configuration includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerTools.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Debug {

   void StorageProfilerTools::printInfo()
   {
      // Analyze the data
      Array reqs;
      Array min;
      Array max;
      StorageProfilerMacro::analyze(reqs, min, max);

      if(QuICCEnv().allowsIO())
      {
         int digits = 1;

         // Create nice looking ouput header
         Tools::Formatter::printNewline(std::cout);
         Tools::Formatter::printLine(std::cout, '%');
         Tools::Formatter::printCentered(std::cout, "Storage profiling information", '%');
         Tools::Formatter::printLine(std::cout, '%');

         // Output the timinigs
         std::stringstream oss;
         int base = 27 + static_cast<int>(min.size() > 0)*10;
         auto namedStorage = StorageProfilerMacro::namePoints();
         std::string memExt;
         MHDFloat memFactor = 0.0;
         Array totMem = Array::Zero(1 + 2*(max.size()>0));
         int i = 0;
         for(auto it = namedStorage.cbegin(); it != namedStorage.cend(); ++it)
         {
            int rem = static_cast<int>(it->first) % StorageProfilerMacro::LVL0;
            if(rem > 0)
            {
               oss << "   ";
               rem = rem % StorageProfilerMacro::LVL1;

               if(rem > 0)
               {
                  oss << "   ";
                  rem = rem % StorageProfilerMacro::LVL2;

                  if(rem > 0)
                  {
                     oss << "   ";
                     rem = rem % StorageProfilerMacro::LVL3;
                  }
               }
            }

            StorageProfilerTools::setUnit(reqs(i), memExt, memFactor);
            oss << it->second << ": " << std::fixed << std::setprecision(digits) << reqs(i)/memFactor;

            if(max.size() != 0)
            {
               oss << " / " << std::fixed << std::setprecision(digits) << max(i)/memFactor << " / " << std::fixed << std::setprecision(digits) << min(i)/memFactor;
            }

            oss << memExt;
            Tools::Formatter::printCentered(std::cout, oss.str(), ' ', base);
            oss.str("");

            // Compute total memory
            if(static_cast<int>(it->first) % StorageProfilerMacro::LVL0 == 0)
            {
               totMem(0) += reqs(i);

               if(totMem.size() > 1)
               {
                  totMem(1) += max(i);
                  totMem(2) += min(i);
               }
            }

            i++;
         }

         Tools::Formatter::printLine(std::cout, '-');
         StorageProfilerTools::setUnit(totMem(0), memExt, memFactor);
         oss << "Total: " << std::fixed << std::setprecision(digits) << totMem(0)/memFactor;
         if(totMem.size() > 1)
         {
            oss << " / " << std::fixed << std::setprecision(digits) << totMem(1)/memFactor << " / " << std::fixed << std::setprecision(digits) << totMem(2)/memFactor;
         }
         oss << memExt;
         Tools::Formatter::printCentered(std::cout, oss.str(), ' ', base);

         Tools::Formatter::printLine(std::cout, '%');
         Tools::Formatter::printNewline(std::cout);
      }
   }

   void StorageProfilerTools::setUnit(MHDFloat value, std::string &ext, MHDFloat &factor)
   {
      if(value < 1024.)
      {
         ext = " B";
         factor = 1.0;
      } else if(value < std::pow(1024.,2))
      {
         ext = " k";
         factor = 1024.;
      } else if(value < std::pow(1024.,3))
      {
         ext = " M";
         factor = 1024.*1024.;
      } else if(value < std::pow(1024.,4))
      {
         ext = " G";
         factor = std::pow(1024.,3);
      }
   }

}
}
