/** 
 * @file TriangularSHmIndexCounter.cpp
 * @brief Source of spherical harmonic index counter with m spectral ordering with Triangular radial truncation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Resolutions/Tools/TriangularSHmIndexCounter.hpp"

// Project includes
//

namespace QuICC {

   TriangularSHmIndexCounter::TriangularSHmIndexCounter(SharedCSimulationResolution spSim, SharedCCoreResolution spCpu)
      : IResolutionIndexCounter(spSim, spCpu)
   {
   }

   TriangularSHmIndexCounter::~TriangularSHmIndexCounter()
   {
   }

   std::vector<int> TriangularSHmIndexCounter::makeVKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const
   {
      std::vector<int> key;

      if(id == Dimensions::Transform::TRA1D)
      {
         key.push_back(i);
         key.push_back(j);
         key.push_back(k);
      } else if(id == Dimensions::Transform::TRA2D)
      {
         key.push_back(j);
         key.push_back(i);
         key.push_back(k);
      } else if(id == Dimensions::Transform::TRA3D)
      {
         key.push_back(k);
         key.push_back(j);
         key.push_back(i);
      }

      return key;
   }

   std::tuple<int,int,int> TriangularSHmIndexCounter::makeKey(const Dimensions::Transform::Id id, const int i, const int j, const int k) const
   {
      std::tuple<int,int,int> key;

      if(id == Dimensions::Transform::TRA1D)
      {
         key = std::make_tuple(i, j, k);
      } else if(id == Dimensions::Transform::TRA2D)
      {
         key = std::make_tuple(j, i, k);
      } else if(id == Dimensions::Transform::TRA3D)
      {
         key = std::make_tuple(k, j, i);
      }

      return key;
   }

   int TriangularSHmIndexCounter::dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId, const MHDFloat) const
   {
      return this->mspSim->dim(simId,spaceId);
   }

   ArrayI TriangularSHmIndexCounter::orderedDimensions(const Dimensions::Space::Id spaceId) const
   {
      ArrayI dims = this->mspSim->dimensions(spaceId);

      return this->orderedDimensions(dims, spaceId);
   }

   ArrayI TriangularSHmIndexCounter::orderedDimensions(const ArrayI& dims, const Dimensions::Space::Id spaceId) const
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
      
      //  Physical space ordering is 3D, 2D, 1D
      } else //if(spaceId == Dimensions::Space::PHYSICAL)
      {
         for(int i = 0; i < dims.size(); ++i)
         {
            oDims(i) = dims(dims.size()-1-i);
         }
      }

      return oDims;
   }

   void TriangularSHmIndexCounter::computeOffsets(std::vector<TriangularSHmIndexCounter::OffsetType>& blocks, std::vector<std::vector<TriangularSHmIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId) const
   {
      this->computeOffsets(blocks, offsets, spaceId, this->mspSim);
   }

   void TriangularSHmIndexCounter::computeOffsets(std::vector<TriangularSHmIndexCounter::OffsetType>& blocks, std::vector<std::vector<TriangularSHmIndexCounter::OffsetType> >& offsets, const Dimensions::Space::Id spaceId, SharedCSimulationResolution spRef) const
   {
      Dimensions::Transform::Id transId;
      Dimensions::Simulation::Id simId;
      
      // Clear the vector of offsets
      offsets.clear();
      std::vector<OffsetType>  offV;

      // In spectral space offset computation, spherical harmonic triangular truncation make it complicated
      if(spaceId == Dimensions::Space::SPECTRAL)
      {
         transId = Dimensions::Transform::TRA1D;
         simId = Dimensions::Simulation::SIM3D;

         // Loop over all local harmonic order m 
         OffsetType offset = 0;
         int m0 = 0;
         for(int iM = 0; iM < this->mspCpu->dim(transId)->dim<Dimensions::Data::DAT3D>(); ++iM)
         {
            int m_ = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT3D>(iM);
            if(m_ < spRef->dim(simId,spaceId))
            {
               // Compute the offset to the local harmonic degree m - 1
               for(int m = m0; m < m_; ++m)
               {  
                  for(int l = m; l < spRef->dim(Dimensions::Simulation::SIM2D,spaceId); ++l)
                  {
                     offset++;
                  }  
               }

               // Compute offset for the local l
               offV.clear();
               for(int iL = 0; iL < this->mspCpu->dim(transId)->dim<Dimensions::Data::DAT2D>(iM); ++iL)
               {
                  int l_ = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT2D>(iL, iM);
                  if(l_ < spRef->dim(Dimensions::Simulation::SIM2D,spaceId))
                  {
                     offV.push_back(offset + l_ - m_);
                  }
               }

               offsets.push_back(offV);

               m0 = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT3D>(iM);

               // 1D blocks
               blocks.push_back(std::min(this->dim(Dimensions::Simulation::SIM1D, spaceId, m_), spRef->dim(Dimensions::Simulation::SIM1D,spaceId)));
            }
         }

      //  Physical space offset computation (regular)
      } else //if(spaceId == Dimensions::Space::PHYSICAL)
      {
         transId = Dimensions::Transform::TRA3D;
         simId = Dimensions::Simulation::SIM1D;

         offV.push_back(0);
         offV.push_back(0);
         for(int i=0; i < this->mspCpu->dim(transId)->dim<Dimensions::Data::DAT3D>(); ++i)
         {
            int i_ = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT3D>(i);
            // Check if value is available in file
            if(i_ < spRef->dim(simId,spaceId))
            {
               // Compute offset for third dimension
               offV.at(0) = i_;

               // Compute offset for second dimension
               offV.at(1) = this->mspCpu->dim(transId)->idx<Dimensions::Data::DAT2D>(0,i);

               offsets.push_back(offV);

               // 1D blocks
               blocks.push_back(std::min(this->dim(Dimensions::Simulation::SIM1D, spaceId, i_), spRef->dim(Dimensions::Simulation::SIM1D,spaceId)));
            }
         }
      }
   }
}
