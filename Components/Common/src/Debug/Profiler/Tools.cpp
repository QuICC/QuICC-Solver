/** 
 * @file Tools.cpp
 * @brief Source of the implementation of tools for the profiling timer
 */

// Configuration includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/Tools.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Typedefs.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Framework/MpiFramework.hpp"

namespace QuICC {

namespace Debug {

namespace Profiler {

   void Tools::writeTimings(std::ostream& stream)
   {
      #ifdef QUICC_PROFILE_PERCORE
         Array ts;
         ProfilerMacro::getTimings(ts);

         // Output communcation structure
         Tools::Formatter::printNewline(stream);
         Tools::Formatter::printLine(stream, '%');
         Tools::Formatter::printCentered(stream, "Communication structure", '%');
         Tools::Formatter::printLine(stream, '%');

         Tools::Formatter::printNewline(stream);
         Tools::Formatter::printCentered(stream, "First data exchange", '-');
         stream << MpiFramework::transformCpus(0).transpose() << std::endl;
         Tools::Formatter::printCentered(stream, "Second data exchange", '-');
         stream << MpiFramework::transformCpus(1).transpose() << std::endl;

         int digits = 3;

         // Create nice looking ouput header
         Tools::Formatter::printNewline(stream);
         Tools::Formatter::printLine(stream, '%');
         Tools::Formatter::printCentered(stream, "Profiling information", '%');
         Tools::Formatter::printLine(stream, '%');

         // Output the timinigs
         std::stringstream oss;
         int base = 27;
         auto namedBreak = ProfilerMacro::namePoints();
         int i = 0;
         for(auto it = namedBreak.begin(); it != namedBreak.end(); ++it)
         {
            oss << it->second << ": " << std::fixed << std::setprecision(digits) << ts(i) << " s";
            Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
            oss.str("");
            ++i;
         }

         Tools::Formatter::printLine(stream, '%');
         Tools::Formatter::printNewline(stream);
      #endif //QUICC_PROFILE_PERCORE
   }

   void Tools::printInfo(std::ostream* pStream)
   {
      if(pStream == nullptr)
      {
         pStream = &std::cout;
      }

      #ifdef QUICC_PROFILE_PERCORE
         Tools::writeTimings(*pStream);
      #endif //QUICC_PROFILE_PERCORE

      // Analyze the data
      Array ts;
      Array min;
      Array max;

      ProfilerMacro::analyze(ts, min, max);

      if(QuICCEnv().allowsIO())
      {
         int digits = 3;

         // Create nice looking ouput header
         QuICC::Tools::Formatter::printNewline(std::cout);
         QuICC::Tools::Formatter::printLine(std::cout, '%');
         QuICC::Tools::Formatter::printCentered(std::cout, "Profiling information", '%');
         QuICC::Tools::Formatter::printLine(std::cout, '%');

         // Output the timinigs
         std::stringstream oss;
         int base = 35 + static_cast<int>(min.size() > 0)*10;
         auto namedBreak = ProfilerMacro::namePoints();
         int i = 0;
         for(auto it = namedBreak.begin(); it != namedBreak.end(); ++it)
         {
            int rem = static_cast<int>(it->first) % LVL0;
            if(rem > 0)
            {
               oss << "   ";
               rem = rem % LVL1;

               if(rem > 0)
               {
                  oss << "   ";
                  rem = rem % LVL2;

                  if(rem > 0)
                  {
                     oss << "   ";
                     rem = rem % LVL3;
                  }
               }
            }

            oss << it->second << ": " << std::scientific << std::setprecision(digits) << ts(i);

            if(max.size() != 0)
            {
               oss << " / " << std::scientific << std::setprecision(digits) << max(i) << " / " << std::scientific << std::setprecision(digits) << min(i);
            }

            oss << " s";
            QuICC::Tools::Formatter::printCentered(std::cout, oss.str(), ' ', base);
            oss.str("");
            i++;
         }

         QuICC::Tools::Formatter::printLine(std::cout, '%');
         QuICC::Tools::Formatter::printNewline(std::cout);
      }
   }

}
}
}
