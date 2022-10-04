/**
 * @file SHm.cpp
 * @brief Source of the tools for schemes with spherical harmonic expansions with m spectral ordering
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Equations/Tools/SHm.hpp"
#include "QuICC/Enums/DimensionTools.hpp"
#include "QuICC/Enums/Dimensions.hpp"

// Project includes
//

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
      for(int m = 0; m < tRes.dim<Dimensions::Data::DAT3D>(); m++)
      {
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(m);
            rTauNs(m) = rTauNs(m)*(res.sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rTauNs(m) = rTauNs(m)*tRes.dim<Dimensions::Data::DAT2D>(m);
         #endif //QUICC_MPISPSOLVE

      }
   }

   void SHm::interpretGalerkinN(ArrayI& rGalerkinNs, const Resolution& res) const
   {
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      for(int m = 0; m < tRes.dim<Dimensions::Data::DAT3D>(); m++)
      {
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(m);
            rGalerkinNs(m) = rGalerkinNs(m)*(res.sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rGalerkinNs(m) = rGalerkinNs(m)*tRes.dim<Dimensions::Data::DAT2D>(m);
         #endif //QUICC_MPISPSOLVE
      }
   }

   void SHm::interpretRhsN(ArrayI&, const Resolution&) const
   {
      // Python setup is sufficient
   }

   void SHm::interpretSystemN(ArrayI& rSystemNs, const Resolution& res) const
   {
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      for(int m = 0; m < tRes.dim<Dimensions::Data::DAT3D>(); m++)
      {
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            int m_ = tRes.idx<Dimensions::Data::DAT3D>(m);
            rSystemNs(m) = rSystemNs(m)*(res.sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rSystemNs(m) = rSystemNs(m)*tRes.dim<Dimensions::Data::DAT2D>(m);
         #endif //QUICC_MPISPSOLVE
      }
   }

}
}
}
