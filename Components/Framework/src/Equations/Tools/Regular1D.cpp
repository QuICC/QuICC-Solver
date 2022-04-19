/**
 * @file Regular1D.cpp
 * @brief Source of the tools for schemes with a single eigen direction
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/Tools/Regular1D.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Equations {

namespace Tools {

   std::vector<MHDFloat> Regular1D::identifyIndexes(const Resolution& res, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      MHDFloat k;
      if(res.sim().ss().dimension() == 3)
      {
         k = static_cast<MHDFloat>(res.cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx));
      } else if(res.sim().ss().dimension() == 2)
      {
         k = static_cast<MHDFloat>(res.cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(matIdx));
      } else
      {
         throw std::logic_error("This eigen direction is not possible for 1D case");
      }

      // Get wave number rescale to box size
      eigs.push_back(res.sim().boxScale(Dimensions::Simulation::SIM2D)*k);

      return eigs;
   }

   int Regular1D::computeNMat(const Resolution& res) const
   {
      int nMat;
      if(res.sim().ss().dimension() == 3)
      {
         nMat = res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      } else if(res.sim().ss().dimension() == 2)
      {
         nMat = res.cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>();
      } else
      {
         throw std::logic_error("This eigen direction is not possible for 1D case");
      }

      return nMat;
   }

   void Regular1D::interpretTauN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void Regular1D::interpretGalerkinN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void Regular1D::interpretRhsN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void Regular1D::interpretSystemN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }
}
}
}
