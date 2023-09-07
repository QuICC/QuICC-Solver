/**
 * @file FullyCoupled.cpp
 * @brief Source of the tools for schemes with no eigen direction
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/Tools/FullyCoupled.hpp"

namespace QuICC {

namespace Equations {

namespace Tools {

   std::vector<MHDFloat> FullyCoupled::identifyIndexes(const Resolution&, const int) const
   {
      std::vector<MHDFloat> eigs;

      return eigs;
   }

   int FullyCoupled::computeNMat(const Resolution&) const
   {
      return 1;
   }

   void FullyCoupled::interpretTauN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void FullyCoupled::interpretGalerkinN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void FullyCoupled::interpretRhsN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void FullyCoupled::interpretSystemN(ArrayI&, const Resolution&, const int) const
   {
      // Python setup is sufficient
   }
}
}
}
