/**
 * @file Tools.cu
 * @brief Defines some useful constants and tools for cuFFT
 */

// System includes
//
#include <cmath>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/Tools.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   const double Tools::STD_DEALIASING = 3.0/2.0;

   const double Tools::COS_DEALIASING = 3.0/2.0;

   const double Tools::MIXED_DEALIASING = 3.0;

   const double Tools::OPTIMIZATION_WIDTH = 0.05;

   bool Tools::optimizeFft(const int size, int& opt)
   {
      // Store the input size
      int factorised;

      // optimisation width (maximum of 5% extension)
      int width = std::max(3, static_cast<int>(Tools::OPTIMIZATION_WIDTH*size));

      bool good = false;
      opt = -1;

      // Loop while optimisation is possible
      while(!good && opt < width)
      {
         // Increment optimisation
         opt++;

         // Get trial size
         factorised = size + opt;

         // Compute factorisation (trial division algorithm)
         while(factorised > 1)
         {
            // Check for optimised factor 2
            if(factorised % 2 == 0)
            {
               factorised /= 2;
               good = true;
            // Check for optimised factor 3
            } else if(factorised % 3 == 0)
            {
               factorised /= 3;
               good = true;
            // Check for optimised factor 5
            } else if(factorised % 5 == 0)
            {
               factorised /= 5;
               good = true;
            // Check for optimised factor 7
            } else if(factorised % 7 == 0)
            {
               factorised /= 7;
               good = true;
            } else
            {
               good = false;
               break;
            }
         }
      }

      return good;
   }

}
}
}
}
}
