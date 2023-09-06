/**
 * @file Regular3D.cpp
 * @brief Source of the tools for schemes with three eigen direction
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Equations/Tools/Regular3D.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

namespace Tools {

   std::vector<MHDFloat> Regular3D::identifyIndexes(const Resolution&, const int) const
   {
      throw std::logic_error("Not yet implemented!");
      std::vector<MHDFloat> eigs;

      // Fill eigs somehow

      return eigs;
   }

   int Regular3D::computeNMat(const Resolution& res) const
   {
      int nMat = 0;

      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += tRes.dim<Dimensions::Data::DAT2D>(i)*tRes.dim<Dimensions::Data::DATB1D>(i);
      }

      return nMat;
   }

   void Regular3D::interpretTauN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void Regular3D::interpretGalerkinN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void Regular3D::interpretRhsN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void Regular3D::interpretSystemN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

}
}
}
