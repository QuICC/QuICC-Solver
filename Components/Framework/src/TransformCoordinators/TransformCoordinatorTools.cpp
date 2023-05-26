/**
 * @file TransformCoordinatorTools.cpp
 * @brief Source of the requirement tools to work with variables and equations
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/TransformCoordinators/TransformCoordinatorTools.hpp"
#include "QuICC/Timers/StageTimer.hpp"

namespace QuICC {

namespace Transform {

   void TransformCoordinatorTools::init(TransformCoordinatorType& rCoord, SharedIForwardGrouper spFwdGrouper, SharedIBackwardGrouper spBwdGrouper, const std::vector<ArrayI>& packs, SharedResolution spRes, const std::map<std::size_t,NonDimensional::SharedINumber>& runOptions)
   {
      StageTimer stage;
      stage.start("initializing transforms");

      // Initialise the transform coordinator
      rCoord.setScheme(spRes->sim().spSpatialScheme());
      for(int i = 0; i < spRes->sim().ss().dimension(); i++)
      {
         Dimensions::Transform::Id id = static_cast<Dimensions::Transform::Id>(i);
         rCoord.addTransform(id, spRes->sim().ss().createTransform(id, spRes->spTransformSetup(id)));
      }

      // Get the list of required options
      std::set<std::size_t>  requests;
      rCoord.requiredOptions(requests);

      // Storage for the options' values
      std::map<std::size_t,NonDimensional::SharedINumber> options;

      // Extract required options
      std::map<std::size_t,NonDimensional::SharedINumber>::const_iterator it;
      for(it = runOptions.begin(); it != runOptions.end(); ++it)
      {
         if(requests.count(it->first) > 0)
         {
            options.insert(std::make_pair(it->first, it->second));
         }
      }

      // Make sure everything was present
      if(requests.size() != options.size())
      {
         throw std::logic_error("TransformCoordinatorTools: requested options have not been set!");
      }

      // Transfer options to transforms
      rCoord.setOptions(options);

      stage.done();
      stage.start("initializing communicators");

      // Initialise the communicator
      rCoord.initCommunicator(spRes);

      stage.done();
      stage.start("initializing converters");

      // Initialise the converters
      rCoord.communicator().initConverter(spRes, packs, spFwdGrouper->split);

      stage.done();
   }

   void TransformCoordinatorTools::addSet(std::vector<ArrayI>& packs, const std::set<int>& packSet)
   {
      packs.push_back(ArrayI(packSet.size()));
      ArrayI& pack = packs.back();
      int pi = 0;
      for(auto k : packSet)
      {
         pack(pi) = k;
         pi++;
      }
   }

   void TransformCoordinatorTools::computePacks(std::vector<ArrayI>& packs, SharedIForwardGrouper spFwdGrouper, SharedIBackwardGrouper spBwdGrouper, const std::map<int, std::vector<Transform::TransformTree> >& forwardTrees, const std::map<int, std::vector<Transform::TransformTree> >& backwardTrees, const std::set<int>& keys,  SharedResolution spRes)
   {
      // Compute pack sizes
      std::set<int>  fwdSet;
      std::set<int>  bwdSet;
      // Get the buffer pack sizes for first dimension
      for(auto j: keys)
      {
         fillSet<stage_1D>(fwdSet, spFwdGrouper, forwardTrees.at(j));

         fillSet<stage_1D>(bwdSet, spBwdGrouper, backwardTrees.at(j));
      }

      // Add sets as packs
      addSet(packs, fwdSet);
      addSet(packs, bwdSet);

      // second dimension
      if(spRes->sim().ss().dimension() == 3)
      {
         fwdSet.clear();
         bwdSet.clear();

         // Get the buffer pack sizes for second dimension
         for(auto j: keys)
         {
            // Fill pack sets
            fillSet<stage_2D>(fwdSet, spFwdGrouper, forwardTrees.at(j));
            fillSet<stage_2D>(bwdSet, spBwdGrouper, backwardTrees.at(j));
         }

         // Add sets as packs
         addSet(packs, fwdSet);
         addSet(packs, bwdSet);
      }
   }
} // Transform
} // QuICC
