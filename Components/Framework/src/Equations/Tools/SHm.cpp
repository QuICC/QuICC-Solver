/**
 * @file SHm.cpp
 * @brief Source of the tools for schemes with spherical harmonic expansions with m spectral ordering
 */

// System includes
//

// Project includes
//
#include "QuICC/Equations/Tools/SHm.hpp"
#include "QuICC/Enums/DimensionTools.hpp"
#include "QuICC/Enums/Dimensions.hpp"

namespace QuICC {

namespace Equations {

namespace Tools {

   std::vector<MHDFloat> SHm::identifyIndexes(const Resolution& res, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      eigs.push_back(static_cast<MHDFloat>(res.cpu()->dim(Dimensions::Transform::SPECTRAL)->idx<Dimensions::Data::DAT3D>(matIdx)));

      return eigs;
   }

   int SHm::computeNMat(const Resolution& res) const
   {
      return res.cpu()->dim(Dimensions::Transform::SPECTRAL)->dim<Dimensions::Data::DAT3D>();
   }

   void SHm::interpretTauN(ArrayI& rTauNs, const Resolution& res) const
   {
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
      {
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            int m = tRes.idx<Dimensions::Data::DAT3D>(k);
            rTauNs(k) = rTauNs(k)*(res.sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m);
         #else
            assert( rTauNs(k) == tRes.dim<Dimensions::Data::DATB1D>(0,k) );
            for(int j = 1; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               rTauNs(k) += tRes.dim<Dimensions::Data::DATB1D>(j,k);
            }
         #endif //QUICC_MPISPSOLVE

      }
   }

   void SHm::interpretGalerkinN(ArrayI& rGalerkinNs, const Resolution& res) const
   {
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
      {
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            int m = tRes.idx<Dimensions::Data::DAT3D>(k);
            rGalerkinNs(k) = rGalerkinNs(k)*(res.sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m);
         #else
            // Get change of truncation due to Galerkin basis
            int gDelta = tRes.dim<Dimensions::Data::DATB1D>(0,k) - rGalerkinNs(k);
            assert( gDelta >= 0 );

            for(int j = 1; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               rGalerkinNs(k) += tRes.dim<Dimensions::Data::DATB1D>(j,k) - gDelta;
            }
         #endif //QUICC_MPISPSOLVE
      }
   }

   void SHm::interpretRhsN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void SHm::interpretSystemN(ArrayI& rSystemNs, const Resolution& res, const int nFields) const
   {
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
      {
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            int m = tRes.idx<Dimensions::Data::DAT3D>(k);
            rSystemNs(k) = rSystemNs(k)*(res.sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m);
         #else
            // Get change of truncation due to Galerkin basis
            int gDelta = nFields*tRes.dim<Dimensions::Data::DATB1D>(0,k) - rSystemNs(k);
            assert( gDelta >= 0 );

            for(int j = 1; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
            {
               rSystemNs(k) += nFields*tRes.dim<Dimensions::Data::DATB1D>(j,k) - gDelta;
            }
         #endif //QUICC_MPISPSOLVE
      }
   }

}
}
}
