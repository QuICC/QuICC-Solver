/**
 * @file Tools.cpp
 * @brief Defines some useful constants and tools for FFT
 */

// System includes
//
#include <cmath>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Tools.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

   int Tools::dealiasFft(const int size)
   {
      return std::ceil(Tools::Backend::STD_DEALIASING*static_cast<MHDFloat>(size));
   }

   int Tools::dealiasComplexFft(const int size)
   {
      return Tools::dealiasFft(size);
   }

   int Tools::dealiasMixedFft(const int size)
   {
      return std::ceil(Tools::Backend::MIXED_DEALIASING*static_cast<MHDFloat>(size));
   }

   int Tools::dealiasCosFft(const int size)
   {
      return std::ceil(Tools::Backend::COS_DEALIASING*static_cast<MHDFloat>(size));
   }

   int Tools::optimizeFft(const int size)
   {
      int opt;
      bool good = Tools::Backend::optimizeFft(size, opt);

      // Size is not good and optimisation failed
      if(!good)
      {
         if(QuICCEnv().allowsIO())
         {
            // Create nice looking warnig message
            QuICC::Tools::Formatter::printLine(std::cout, '%');
            QuICC::Tools::Formatter::printCentered(std::cout, "WARNING: FFT size optimization failed", '%');
            QuICC::Tools::Formatter::printCentered(std::cout, "Selected FFT size might be slow!", '%');
            QuICC::Tools::Formatter::printLine(std::cout, '%');
            QuICC::Tools::Formatter::printNewline(std::cout);
         }

         return size;

      // Size was not good but optimisation worked
      } else
      {
         // Create nice looking warning message
         if(opt > 0)
         {
            if(QuICCEnv().allowsIO())
            {
               QuICC::Tools::Formatter::printLine(std::cout, '%');
               std::stringstream oss;
               if(opt > 2)
               {
                  oss << "WARNING: ";
               }
               oss << "Extended FFT size (+" << opt << ")!";
               QuICC::Tools::Formatter::printCentered(std::cout, oss.str(), '%');
               QuICC::Tools::Formatter::printLine(std::cout, '%');
               QuICC::Tools::Formatter::printNewline(std::cout);
            }
         }

         return size + opt;
      }
   }

}
}
}
