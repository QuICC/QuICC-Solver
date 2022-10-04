/** 
 * @file IRegularSHmlBuilder.cpp
 * @brief Source of the Regular basis + Spherical Harmonics scheme implementation with m spectral ordering
 */

// System includes
//
#include <cassert>
#include <set>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/SpatialScheme/3D/IRegularSHmlBuilder.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/Tools/Regular.hpp"
#include "QuICC/SpatialScheme/Tools/SH.hpp"
#include "QuICC/Resolutions/Tools/RegularSHmlIndexCounter.hpp"

namespace QuICC {

namespace SpatialScheme {

   IRegularSHmlBuilder::IRegularSHmlBuilder(const ArrayI& dim, const GridPurpose::Id purpose)
      : IBuilder(dim.size(), purpose), mI(dim(0)), mL(dim(1)), mM(dim(2))
   {
      assert(dim.size() == 3);
   }

   IRegularSHmlBuilder::~IRegularSHmlBuilder()
   {
   }

   ArrayI IRegularSHmlBuilder::resolution() const
   {
      ArrayI space(this->dims());
      space << this->mI, this->mL, this->mM;

      return space;
   }

   void IRegularSHmlBuilder::addIndexCounter(SharedResolution spRes)
   {
      auto spCounter = std::make_shared<RegularSHmlIndexCounter>(spRes->spSim(), spRes->spCpu());

      spRes->setIndexCounter(spCounter);
   }

   int IRegularSHmlBuilder::fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Splitting::Locations::Id flag)
   {
      // Safety assertions for default values
      assert( (id.size() == 0) || (bins.size() > 0) );
      assert( id.size() == bins.size() );
      assert( n0.size() == nN.size() );
      assert( (bins.size() == 0) || (flag != Splitting::Locations::NONE) );

      // Make sure we start with empty indexes
      fwd1D.clear();
      bwd1D.clear();
      idx2D.clear();

      // Multimap for the modes
      std::multimap<int,int> modes;

      // No splitting
      if(flag == Splitting::Locations::NONE)
      {
         this->splitSerial(modes, transId);
      }
      // Splitting is on first transform
      else if(flag == Splitting::Locations::FIRST)
      {
         this->splitSingle1D(modes, n0, nN, transId);
      }
      // Splitting is on second transform
      else if(flag == Splitting::Locations::SECOND)
      {
         this->splitSingle2D(modes, id, bins, n0, nN, transId);
      }
      // Splitting is on both transforms
      else if(flag == Splitting::Locations::BOTH)
      {
         this->splitTubular(modes, id, bins, n0, nN, transId);
      }

      // Fill indexes for 2D and 3D
      Tools::Regular::fillIndexes2D3D(idx2D, idx3D, modes);

      // Fill indexes for 1D
      if(transId == Dimensions::Transform::TRA1D || transId == Dimensions::Transform::TRA3D || transId == Dimensions::Transform::SPECTRAL)
      {
         Tools::Regular::fillIndexes1D(fwd1D, bwd1D, idx3D, this->dim(transId, Dimensions::Data::DATF1D), this->dim(transId, Dimensions::Data::DATB1D));
      }
      else if(transId == Dimensions::Transform::TRA2D)
      {
         Tools::SH::fillIndexes1D(fwd1D, bwd1D, idx3D, this->dim(transId, Dimensions::Data::DATF1D), this->dim(transId, Dimensions::Data::DATB1D));
      }

      // Set status (0 for success, 1 for failure)
      int status = 0;
      if(modes.size() == 0 && nN.size() > 0)
      {
         status = 1;
      }

      return status;
   }

   int IRegularSHmlBuilder::splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag)
   {
      // Splittable size for first transform splitting
      if(flag == Splitting::Locations::FIRST)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);
            int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

            return Tools::SH::nHarmonics(nL, nM);
         }
         // Get total size for second transform
         else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);
         }
         // Get total size for third transform
         else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }
         // Get total size for spectral space
         else if(transId == Dimensions::Transform::SPECTRAL)
         {
            int nL = this->dim(transId, Dimensions::Data::DAT2D);
            int nM = this->dim(transId, Dimensions::Data::DAT3D);

            return Tools::SH::nHarmonics(nL, nM);
         }
      }
      // Splittable size for second transform splitting
      else if(flag == Splitting::Locations::SECOND)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);
         }
         // Get total size for second transform
         else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }
         // Get total size for third transform
         else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);
         }
         // Get total size for spectral space
         else if(transId == Dimensions::Transform::SPECTRAL)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }
      }
      // Splittable size for both transforms splitting
      else if(flag == Splitting::Locations::BOTH)
      {
         // Get total size for first transform
         if(transId == Dimensions::Transform::TRA1D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }
         // Get total size for second transform
         else if(transId == Dimensions::Transform::TRA2D)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }
         // Get total size for third transform
         else if(transId == Dimensions::Transform::TRA3D)
         {
            return this->dim(transId, Dimensions::Data::DAT2D);
         }
         // Get total size for spectral space
         else if(transId == Dimensions::Transform::SPECTRAL)
         {
            return this->dim(transId, Dimensions::Data::DAT3D);
         }
      }
      
      // If none of the previous conditions were right
      throw std::logic_error("Tried to split in a unknown dimension for regular spherical harmonics case");

      return -1;
   }

   void IRegularSHmlBuilder::splitSerial(std::multimap<int,int>& modes, const Dimensions::Transform::Id transId)
   {
      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT3D);
         int nM = this->dim(transId, Dimensions::Data::DAT2D);

         // Get full list of harmonics mapped by harmonic degree l
         Tools::SH::buildLMap(modes, nL, nM);
      }
      // Create index list for spectral space
      else if(transId == Dimensions::Transform::SPECTRAL)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT2D);
         int nM = this->dim(transId, Dimensions::Data::DAT3D);

         // Get full list of harmonics mapped by harmonic order m
         Tools::SH::buildMMap(modes, nL, nM);
      }
      // Create index list for second and third transform
      else
      {
         int k0 = 0;
         int kN = this->dim(transId, Dimensions::Data::DAT3D);
         ArrayI j0 = ArrayI::Zero(kN);
         ArrayI jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
         int c0 = 0;
         int cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
   }

   void IRegularSHmlBuilder::splitSingle1D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT3D);
         int nM = this->dim(transId, Dimensions::Data::DAT2D);

         // Get restricted sorted list of harmonics
         Tools::SH::buildLHSortedMap(modes, nL, nM, n0(0), nN(0));
      }
      // Create index list for second transform
      else if(transId == Dimensions::Transform::TRA2D)
      {
         int k0 = 0;
         int kN = this->dim(transId, Dimensions::Data::DAT3D);
         ArrayI j0 = ArrayI::Constant(kN, n0(0));
         ArrayI jN = ArrayI::Constant(kN, nN(0));
         int c0 = 0;
         int cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
      // Create index list for third transform
      else if(transId == Dimensions::Transform::TRA3D)
      {
         int k0 = n0(0);
         int kN = nN(0);
         ArrayI j0 = ArrayI::Zero(kN);
         ArrayI jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
         int c0 = 0;
         int cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
      // Create index list for spectral space
      else if(transId == Dimensions::Transform::SPECTRAL)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT2D);
         int nM = this->dim(transId, Dimensions::Data::DAT3D);

         // Get restricted sorted list of harmonics
         Tools::SH::buildMHSortedMap(modes, nL, nM, n0(0), nN(0));
      }
   }

   void IRegularSHmlBuilder::splitSingle2D(std::multimap<int,int>& modes, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT3D);
         int nM = this->dim(transId, Dimensions::Data::DAT2D);

         // Get restricted list of harmonics
         std::vector<int> binOrders;
         Tools::SH::binMLLoad(binOrders, nL, nM, id(0), bins(0));

         // Fill with correct modes
         for(auto vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int l = *vIt; l < this->dim(transId, Dimensions::Data::DAT2D); l++)
            {
               modes.insert(std::make_pair(l, *vIt));
            }
         }
      }
      // Create index list for second transform
      else if(transId == Dimensions::Transform::TRA2D)
      {
         int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);
         int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

         // Get restricted list of harmonics
         std::vector<int> binOrders;
         Tools::SH::binMLLoad(binOrders, nL, nM, id(0), bins(0));

         // Fill with correct modes
         for(auto vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int r = 0; r < this->dim(transId, Dimensions::Data::DAT2D); r++)
            {
               modes.insert(std::make_pair(*vIt, r));
            }
         }
      }
      // Create index list for third transform
      else if(transId == Dimensions::Transform::TRA3D)
      {
         int k0 = 0;
         int kN = this->dim(transId, Dimensions::Data::DAT3D);
         ArrayI j0 = ArrayI::Zero(kN);
         ArrayI jN = ArrayI::Constant(kN, this->dim(transId, Dimensions::Data::DAT2D));
         int c0 = n0(0);
         int cN = n0(0) + nN(0);

         Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
      // Create index list for spectral space
      else if(transId == Dimensions::Transform::SPECTRAL)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT2D);
         int nM = this->dim(transId, Dimensions::Data::DAT3D);

         // Get restricted list of harmonics
         std::vector<int> binOrders;
         Tools::SH::binMLLoad(binOrders, nL, nM, id(0), bins(0));

         // Fill with correct modes
         for(auto vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int l = *vIt; l < this->dim(transId, Dimensions::Data::DAT3D); l++)
            {
               modes.insert(std::make_pair(*vIt, l));
            }
         }
      }
   }

   void IRegularSHmlBuilder::splitTubular(std::multimap<int,int>& modes, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId)
   {
      // Create index list for first transform
      if(transId == Dimensions::Transform::TRA1D)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT3D);
         int nM = this->dim(transId, Dimensions::Data::DAT2D);

         // Get restricted list of harmonic orders
         std::vector<int> binOrders;
         Tools::SH::binMLLoad(binOrders, nL, nM, id(1), bins(1));

         // Fill with correct modes
         std::multimap<int,int> tmp;
         for(auto vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int l = *vIt; l < this->dim(transId, Dimensions::Data::DAT2D); l++)
            {
               tmp.insert(std::make_pair(l, *vIt));
            }
         }

         // Get count of balanced harmonic distribution
         ArrayI binLoad = ArrayI::Zero(bins(0));
         for(int i = 0; i < static_cast<int>(tmp.size()); i++)
         {
            binLoad(i % bins(0)) += 1;
         }
         int curLoad = binLoad(id(0));

         std::multimap<int,int> incomplete;
         std::multimap<int,int>::iterator mapIt;
         std::pair<std::multimap<int,int>::iterator, std::multimap<int,int>::iterator> mapRange;
         int curId = 0;
         int counter = 0;
         int l = nL-1; 
         size_t maxLoop = tmp.size() + 5*bins(0);
         for(size_t loop = 0; loop < maxLoop; ++loop)
         {
            if(binLoad(curId) > 0)
            {
               // Extract modes from full list
               if(incomplete.size() == 0 || tmp.count(l) > static_cast<size_t>(incomplete.rbegin()->first))
               {
                  // Get range of m for current l
                  mapRange = tmp.equal_range(l);
                  --l;
               }
               // Get modes from incomplete storage
               else
               {
                  // Get range for largest incomplete set
                  mapRange = tmp.equal_range(incomplete.rbegin()->second);

                  // Delete it from incomplete map
                  mapIt = incomplete.end();
                  --mapIt;
                  incomplete.erase(mapIt);
               }

               // Get number of harmonic orders in range
               int lCount = std::distance(mapRange.first, mapRange.second);

               // Extract all m
               if(lCount <= binLoad(curId))
               {
                  // Store the modes
                  if(curId == id(0))
                  {
                     modes.insert(mapRange.first, mapRange.second);
                  }

                  // Delete used entries
                  tmp.erase(mapRange.first, mapRange.second);

                  // Substract used modes
                  binLoad(curId) -= lCount;
               }
               // Extract partial m
               else
               {
                  // Add information about remaining modes in incomplete list
                  incomplete.insert(std::make_pair(lCount - binLoad(curId), mapRange.first->first));

                  // Create proper range
                  mapRange.second = mapRange.first;
                  std::advance(mapRange.second, binLoad(curId));

                  // Store the modes
                  if(curId == id(0))
                  {
                     modes.insert(mapRange.first, mapRange.second);
                  }

                  // Delete used entries
                  tmp.erase(mapRange.first, mapRange.second);

                  // bin is full
                  binLoad(curId) = 0;
               }
            }

            // Shortcut out of loop
            if(binLoad(id(0)) == 0 || tmp.size() == 0)
            {
               break;
            }

            // Update loop counter
            ++counter;

            // Fill by alternating directions
            curId = counter % bins(0);
            if(counter/bins(0) % 2 == 1)
            {
               curId = bins(0) - curId - 1;
            }

            assert(l >= 0);
         }

         assert(binLoad(id(0)) == 0);
         assert(modes.size() == static_cast<size_t>(curLoad));
      }
      // Create index list for second transform
      else if(transId == Dimensions::Transform::TRA2D)
      {
         int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);
         int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

         // Get restricted list of harmonics
         std::vector<int> binOrders;
         Tools::SH::binMLLoad(binOrders, nL, nM, id(1), bins(1));

         // Fill with correct modes
         for(auto vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int r = 0; r < nN(1); r++)
            {
               modes.insert(std::make_pair(*vIt, n0(1) + r));
            }
         }
      }
      // Create index list for third transform
      else if(transId == Dimensions::Transform::TRA3D)
      {
         int k0 = n0(0);
         int kN = nN(0);
         ArrayI j0 = n0.tail(kN);
         ArrayI jN = nN.tail(kN);
         int c0 = 0;
         int cN = this->dim(transId, Dimensions::Data::DAT2D)*this->dim(transId, Dimensions::Data::DAT3D);

         Tools::Regular::buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
      // Create index list for spectral space
      else if(transId == Dimensions::Transform::SPECTRAL)
      {
         int nL = this->dim(transId, Dimensions::Data::DAT2D);
         int nM = this->dim(transId, Dimensions::Data::DAT3D);

         // Get restricted list of harmonic orders
         std::vector<int> binOrders;
         Tools::SH::binMLLoad(binOrders, nL, nM, id(1), bins(1));

         // Fill with correct modes
         std::multimap<int,int> tmp;
         for(auto vIt = binOrders.begin(); vIt != binOrders.end(); ++vIt)
         {
            for(int l = *vIt; l < this->dim(transId, Dimensions::Data::DAT3D); l++)
            {
               tmp.insert(std::make_pair(*vIt, l));
            }
         }

         int i = 0;
         for(auto mapIt = tmp.begin(); mapIt != tmp.end(); ++mapIt)
         {
            // add mode
            if(i % bins(0) == id(0))
            {
               modes.insert(*mapIt);
            }

            ++i;
         }
      }
   }

   void IRegularSHmlBuilder::buildMLMap(std::multimap<int,int>& harmonics, const int id, const int bins)
   {
      // Assert that more than one bin is available
      assert(bins > 1);

      // Make sure the list of harmonics is empty
      harmonics.clear();

      int nL = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT3D);
      int nM = this->dim(Dimensions::Transform::TRA1D, Dimensions::Data::DAT2D);

      // Reset the loads
      this->resetLoad(); // ICosts

      // Initialise the loads
      Tools::SH::initMLLoad(this->mLoadList, this->mLoad, this->mOptimalLoad, nL, nM, bins);

      // Combine load into pairs
      Tools::SH::combineMPairs(this->mLoadList, this->mRegularLoad, nM, bins);

      // Fill bins with pairs
      Tools::SH::fillMPairsBins(this->mLoadList, this->mRegularLoad, this->mLoad, bins);

      // Update Load sums
      this->updateLoad(bins); // ICosts

      // Fill bins with remaining modes
      Tools::SH::fillMRestBins(this->mLoadList, this->mRegularLoad, this->mLoad, this->mOptimalLoad, bins);

      // Update Load sums
      this->updateLoad(bins); // ICosts

      // Convert loads to harmonic orders
      Tools::SH::convertLoadToOrders(this->mRegularLoad, nL);

      // Extract correct bin
      auto mapRange = this->mRegularLoad.equal_range(id);

      // Loop over all orders
      for(auto mapIt = mapRange.first; mapIt != mapRange.second; mapIt++)
      {
         for(int l = mapIt->second; l < nL; l++)
         {
            harmonics.insert(std::make_pair(mapIt->second, l));
         }
      }
   }
}
}
