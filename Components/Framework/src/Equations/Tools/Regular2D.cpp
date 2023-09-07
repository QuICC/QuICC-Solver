/**
 * @file Regular2D.cpp
 * @brief Source of the tools for schemes with two eigen direction
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/Tools/Regular2D.hpp"
#include "QuICC/Enums/Dimensions.hpp"

namespace QuICC {

namespace Equations {

namespace Tools {

   std::vector<MHDFloat> Regular2D::identifyIndexes(const Resolution& res, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      // Get mode indexes
      ArrayI mode = res.cpu()->dim(Dimensions::Transform::SPECTRAL)->mode(matIdx);

      int sN = res.sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);

      // k2D_
      if(mode(2) < sN/2 + (sN % 2))
      {
         eigs.push_back(res.sim().boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(2)));
      } else
      {
         eigs.push_back(res.sim().boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(2) - sN));
      }

      // k3D_
      eigs.push_back(res.sim().boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(3)));

      return eigs;
   }

   int Regular2D::computeNMat(const Resolution& res) const
   {
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      int nMat = 0;
      for(int i = 0; i < tRes.dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += tRes.dim<Dimensions::Data::DAT2D>(i);
      }

      return nMat;
   }

   void Regular2D::interpretTauN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void Regular2D::interpretGalerkinN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void Regular2D::interpretRhsN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void Regular2D::interpretSystemN(ArrayI&, const Resolution&, const int) const
   {
      // Python setup is sufficient
   }

}
}
}
