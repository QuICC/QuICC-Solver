/** 
 * @file RegularIndexCounter.cpp
 * @brief Source of regular index counter
 */

// System includes
//

// Project includes
//
#include "QuICC/Resolutions/Tools/RegularIndexCounter.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

   RegularIndexCounter::RegularIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu)
      : IResolutionIndexCounter(spSim, spCpu)
   {
   }

   ArrayI RegularIndexCounter::orderedDimensions(const Dimensions::Space::Id spaceId) const
   {
      ArrayI dims = this->mspSim->dimensions(spaceId);

      return this->orderedDimensions(dims, spaceId);
   }

   ArrayI RegularIndexCounter::orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const
   {
      // Storage for the ordered dimensions
      ArrayI oDims(dims.size());

      // Spectral and transform space ordering is 1D, 3D, 2D
      if(spaceId == Dimensions::Space::SPECTRAL || spaceId == Dimensions::Space::TRANSFORM)
      {
         oDims(0) = dims(0);
         for(int i = 1; i < dims.size(); ++i)
         {
            oDims(i) = dims(dims.size()-i);
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

   void RegularIndexCounter::computeOffsets(std::vector<std::vector<RegularIndexCounter::OffsetType>>& blocks, std::vector<std::vector<RegularIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const
   {
      this->computeOffsets(blocks, offsets, spaceId, this->mspSim);
   }

   void RegularIndexCounter::computeOffsets(std::vector<std::vector<RegularIndexCounter::OffsetType>>& blocks, std::vector<std::vector<RegularIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const
   {
      Dimensions::Transform::Id transId;
      Dimensions::Simulation::Id simId;

      // Select transform dimension depending on dimension space
      if(spaceId == Dimensions::Space::SPECTRAL || spaceId == Dimensions::Space::TRANSFORM)
      {
         transId = Dimensions::Transform::SPECTRAL;
         simId = Dimensions::Simulation::SIM2D;

      }
      else //if(spaceId == Dimensions::Space::PHYSICAL)
      {
         transId = static_cast<Dimensions::Transform::Id>(this->mspSim->ss().dimension()-1);
         simId = Dimensions::Simulation::SIM1D;
      }

      // Clear the vector of offsets
      offsets.clear();
      std::vector<OffsetType>  offV;

      auto&& tRes = *this->mspCpu->dim(transId);

      if(this->mspSim->ss().dimension() == 3)
      {
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
               std::vector<RegularIndexCounter::OffsetType> blk;
               blk.push_back(std::min(this->dim(Dimensions::Simulation::SIM1D, spaceId, i_), spRef->dim(Dimensions::Simulation::SIM1D,spaceId)));
               blocks.push_back(blk);
            }
         }
      }
      else if(this->mspSim->ss().dimension() == 2)
      {
         offV.push_back(0);
         offV.push_back(0);
         offV.push_back(0);
         for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(); ++j)
         {
            int j_ = tRes.idx<Dimensions::Data::DAT2D>(j);
            // Check if value is available in file
            if(j_ < spRef->dim(simId,spaceId))
            {
               // Compute offset for third dimension
               offV.at(0) = tRes.idx<Dimensions::Data::DAT2D>(0);

               // Store 2D index
               offV.at(1) = j;

               // Store (fake) 3D index
               offV.at(2) = 0;

               offsets.push_back(offV);

               // 1D blocks
               std::vector<RegularIndexCounter::OffsetType> blk;
               blk.push_back(std::min(this->dim(Dimensions::Simulation::SIM1D, spaceId, j_), spRef->dim(Dimensions::Simulation::SIM1D,spaceId)));
               blocks.push_back(blk);
            }
         }
      }
      else if(this->mspSim->ss().dimension() == 1)
      {
         // Compute offset for third dimension
         offV.at(0) = 0;

         // Store 2D index
         offV.at(1) = 0;

         // Store (fake) 3D index
         offV.at(2) = 0;

         offsets.push_back(offV);

         // 1D blocks
         std::vector<RegularIndexCounter::OffsetType> blk;
         blk.push_back(std::min(this->dim(Dimensions::Simulation::SIM1D, spaceId, 0), spRef->dim(Dimensions::Simulation::SIM1D,spaceId)));
         blocks.push_back(blk);
      }
   }
} // QuICC
