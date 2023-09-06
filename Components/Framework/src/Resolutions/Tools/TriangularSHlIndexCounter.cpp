/**
 * @file TriangularSHlIndexCounter.cpp
 * @brief Source of spherical harmonic index counter with l spectral ordering with triangular radial truncation
 */

// System includes
//

// Project includes
//
#include "QuICC/Resolutions/Tools/TriangularSHlIndexCounter.hpp"
#include "QuICC/SpatialScheme/Tools/TriangularSH.hpp"

namespace QuICC {

   TriangularSHlIndexCounter::TriangularSHlIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu)
      : IResolutionIndexCounter(spSim, spCpu)
   {
   }

   int TriangularSHlIndexCounter::dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId, const MHDFloat l) const
   {
      if(simId == Dimensions::Simulation::SIM1D && spaceId != Dimensions::Space::PHYSICAL)
      {
         SpatialScheme::Tools::TriangularSH t;
         return t.truncationBwd(this->mspSim->dim(simId, spaceId), static_cast<int>(l));
      }
      else
      {
         return this->mspSim->dim(simId, spaceId);
      }
   }

   ArrayI TriangularSHlIndexCounter::orderedDimensions(const Dimensions::Space::Id spaceId) const
   {
      ArrayI dims = this->mspSim->dimensions(spaceId);

      return this->orderedDimensions(dims, spaceId);
   }

   ArrayI TriangularSHlIndexCounter::orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const
   {
      // Storage for the ordered dimensions
      ArrayI oDims(dims.size());

      // In spectral space reduce dimensions to 1D, NH (=number of harmonics)
      if(spaceId == Dimensions::Space::SPECTRAL)
      {
         assert(dims.size() == 3);

         oDims.resize(2);
         oDims(0) = dims(0);
         oDims(1) = 0;

         for(int l = 0; l < dims(1); ++l)
         {
            for(int m = 0; m < std::min(l+1,dims(2)); ++m)
            {
               oDims(1)++;
            }
         }

      }
      //  Physical space ordering is 3D, 2D, 1D
      else //if(spaceId == Dimensions::Space::PHYSICAL)
      {
         for(int i = 0; i < dims.size(); ++i)
         {
            oDims(i) = dims(dims.size()-1-i);
         }
      }

      return oDims;
   }

   void TriangularSHlIndexCounter::computeOffsets(std::vector<TriangularSHlIndexCounter::OffsetType>& blocks, std::vector<std::vector<TriangularSHlIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const
   {
      this->computeOffsets(blocks, offsets, spaceId, this->mspSim);
   }

   void TriangularSHlIndexCounter::computeOffsets(std::vector<TriangularSHlIndexCounter::OffsetType>& blocks, std::vector<std::vector<TriangularSHlIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const
   {
      Dimensions::Simulation::Id simId;

      // Clear the vector of offsets
      offsets.clear();
      std::vector<OffsetType>  offV;

      // In spectral space offset computation, spherical harmonic triangular truncation make it complicated
      if(spaceId == Dimensions::Space::SPECTRAL)
      {
         const auto& tRes = *this->mspCpu->dim(Dimensions::Transform::SPECTRAL);
         simId = Dimensions::Simulation::SIM2D;

         // Loop over all local harmonic degrees l
         OffsetType offset = 0;
         int l0 = 0;
         for(int iL = 0; iL < tRes.dim<Dimensions::Data::DAT3D>(); ++iL)
         {
            int l_ = tRes.idx<Dimensions::Data::DAT3D>(iL);
            if(l_ < spRef->dim(simId,spaceId))
            {
               // Compute the offset to the local harmonic degree l - 1
               for(int l = l0; l < l_; ++l)
               {
                  for(int m = 0; m < std::min(l+1,spRef->dim(Dimensions::Simulation::SIM3D,spaceId)); ++m)
                  {
                     offset++;
                  }
               }

               // Compute offset for the local m
               offV.clear();
               for(int iM = 0; iM < tRes.dim<Dimensions::Data::DAT2D>(iL); ++iM)
               {
                  int m_ = tRes.idx<Dimensions::Data::DAT2D>(iM, iL);
                  if(m_ < spRef->dim(Dimensions::Simulation::SIM3D,spaceId))
                  {
                     offV.push_back(offset + m_);
                  }
               }

               offsets.push_back(offV);

               l0 = tRes.idx<Dimensions::Data::DAT3D>(iL);

               // 1D blocks
               blocks.push_back(std::min(this->dim(Dimensions::Simulation::SIM1D, spaceId, l_), spRef->dim(Dimensions::Simulation::SIM1D,spaceId)));
            }
         }
      }
      // Physical space offset computation (regular)
      else //if(spaceId == Dimensions::Space::PHYSICAL)
      {
         const auto& tRes = *this->mspCpu->dim(Dimensions::Transform::TRA3D);
         simId = Dimensions::Simulation::SIM1D;

         offV.push_back(0);
         offV.push_back(0);
         offV.push_back(0);
         for(int i=0; i < tRes.dim<Dimensions::Data::DAT3D>(); ++i)
         {
            int i_ = tRes.idx<Dimensions::Data::DAT3D>(i);
            // Check if value is available in file
            if(i_ < spRef->dim(simId,spaceId))
            {
               // Compute offset for third dimension
               offV.at(0) = i_;

               // Compute offset for second dimension
               offV.at(1) = tRes.idx<Dimensions::Data::DAT2D>(0,i);

               // Store 3D index
               offV.at(2) = i;

               offsets.push_back(offV);

               // 1D blocks
               blocks.push_back(std::min(this->dim(Dimensions::Simulation::SIM1D, spaceId, i_), spRef->dim(Dimensions::Simulation::SIM1D,spaceId)));
            }
         }
      }
   }
} // QuICC
