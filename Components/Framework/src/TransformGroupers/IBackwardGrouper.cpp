/**
 * @file IBackwardGrouper.cpp
 * @brief Source of the implementation of the equation wise forward transform grouper in xD space
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <set>

// Class include
//
#include "QuICC/TransformGroupers/IBackwardGrouper.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

   IBackwardGrouper::IBackwardGrouper()
      : split(Splitting::Locations::NONE)
   {
   }

   IBackwardGrouper::~IBackwardGrouper()
   {
   }

   int IBackwardGrouper::packs1D(const TransformTree& tree)
   {
      int packs = tree.nEdges(0);

      return packs;
   }

   int IBackwardGrouper::packs2D(const TransformTree& tree)
   {
      int packs = tree.nEdges(1);

      return packs;
   }

   ArrayI IBackwardGrouper::namePacks1D(const std::vector<TransformTree>& projectorTree)
   {
      // Create list of packet sizes
      std::set<int>  list;

      // Loop over all edges
      for(auto treeIt = projectorTree.cbegin(); treeIt != projectorTree.cend(); ++treeIt)
      {
         int counter = this->packs1D(*treeIt);
         list.insert(counter);

         this->mNamedPacks1D.insert(std::make_pair(std::make_pair(treeIt->name(),treeIt->comp<FieldComponents::Spectral::Id>()), counter));
      }

      // Initialise the number of packs
      ArrayI packs(list.size());

      // Set packet sizes
      int i = 0;
      for(auto it = list.cbegin(); it != list.cend(); ++it, ++i)
      {
         packs(i) = *it;
      }

      return packs;
   }

   ArrayI IBackwardGrouper::groupPacks1D(const std::vector<TransformTree>& projectorTree)
   {
      // Initialise the number of packs
      ArrayI packs = this->namePacks1D(projectorTree);

      // Resize packs to single value
      packs.resize(1);
      packs.setConstant(0);

      for(auto it = this->mNamedPacks1D.cbegin(); it != this->mNamedPacks1D.cend(); ++it)
      {
         packs(0) = packs(0) + it->second;
      }

      return packs;
   }

   ArrayI IBackwardGrouper::namePacks2D(const std::vector<TransformTree>& projectorTree)
   {
      // Create list of packet sizes
      std::set<int>  list;

      // Loop over all edges
      for(auto treeIt = projectorTree.cbegin(); treeIt != projectorTree.cend(); ++treeIt)
      {
         int counter = this->packs2D(*treeIt);
         list.insert(counter);

         this->mNamedPacks2D.insert(std::make_pair(std::make_pair(treeIt->name(),treeIt->comp<FieldComponents::Spectral::Id>()), counter));
      }

      // Initialise the number of packs
      ArrayI packs(list.size());

      // Set packet sizes
      int i = 0;
      for(auto it = list.cbegin(); it != list.cend(); ++it, ++i)
      {
         packs(i) = *it;
      }

      return packs;
   }

   ArrayI IBackwardGrouper::groupPacks2D(const std::vector<TransformTree>& projectorTree)
   {
      // Initialise the number of packs
      ArrayI packs = this->namePacks2D(projectorTree);

      // Resize packs to single value
      packs.resize(1);
      packs.setConstant(0);

      for(auto it = this->mNamedPacks2D.cbegin(); it != this->mNamedPacks2D.cend(); ++it)
      {
         packs(0) = packs(0) + it->second;
      }

      return packs;
   }

}
}
