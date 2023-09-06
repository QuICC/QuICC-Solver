/** 
 * @file DoublePeriodicIndexCounter.cpp
 * @brief Source of double periodic index counter
 */

// System includes
//

// Project includes
//
#include "QuICC/Resolutions/Tools/DoublePeriodicIndexCounter.hpp"

namespace QuICC {

   DoublePeriodicIndexCounter::DoublePeriodicIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu)
      : IResolutionIndexCounter(spSim, spCpu)
   {
   }

   ArrayI DoublePeriodicIndexCounter::orderedDimensions(const Dimensions::Space::Id spaceId) const
   {
      ArrayI dims = this->mspSim->dimensions(spaceId);

      return this->orderedDimensions(dims, spaceId);
   }

   ArrayI DoublePeriodicIndexCounter::orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const
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

   void DoublePeriodicIndexCounter::computeOffsets(std::vector<DoublePeriodicIndexCounter::OffsetType>& blocks, std::vector<std::vector<DoublePeriodicIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const
   {
      this->computeOffsets(blocks, offsets, spaceId, this->mspSim);
   }

   void DoublePeriodicIndexCounter::computeOffsets(std::vector<DoublePeriodicIndexCounter::OffsetType>& blocks, std::vector<std::vector<DoublePeriodicIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const
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
         transId = Dimensions::Transform::TRA3D;
         simId = Dimensions::Simulation::SIM1D;
      }

      // Clear the vector of offsets
      offsets.clear();
      std::vector<OffsetType>  offV;

      offV.push_back(0);
      offV.push_back(0);
      offV.push_back(0);

      auto&& tRes = *this->mspCpu->dim(transId);

      // Select transform dimension depending on dimension space
      if(spaceId == Dimensions::Space::SPECTRAL)
      {
         // Get full slowest resolution resolution
         int dat3D = this->mspSim->dim(simId, Dimensions::Space::TRANSFORM);
         int sim3D = this->mspSim->dim(simId,spaceId)/2 + 1;
         int ref3D = spRef->dim(simId,spaceId)/2 + 1;
         int min3D = std::min(sim3D,ref3D);

         for(int i=0; i < tRes.dim<Dimensions::Data::DAT3D>(); ++i)
         {
            int i_ = tRes.idx<Dimensions::Data::DAT3D>(i);
            // Check if value is available in file
            if(i_ < min3D)
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
            else if(dat3D - i_ < min3D)
            {
               // Compute offset for third dimension
               offV.at(0) = i_ - (dat3D - spRef->dim(simId,spaceId));

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
      else //if(spaceId == Dimensions::Space::PHYSICAL || spaceId == Dimensions::Space::TRANSFORM)
      {
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
