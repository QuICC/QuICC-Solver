/**
 * @file Regular2D.cpp
 * @brief Source of the tools for schemes with two eigen direction
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/Tools/Regular2D.hpp"

// Project includes
//

namespace QuICC {

namespace Equations {

namespace Tools {

   std::vector<MHDFloat> Regular2D::identifyIndexes(const Resolution& res, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      // Get mode indexes
      ArrayI mode = res.cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

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
      int nMat = 0;

      for(int i = 0; i < res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
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

   void Regular2D::interpretSystemN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

}
}
}
