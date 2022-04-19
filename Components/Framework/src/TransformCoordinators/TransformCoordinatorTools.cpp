/**
 * @file TransformCoordinatorTools.cpp
 * @brief Source of the requirement tools to work with variables and equations
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/TransformCoordinators/TransformCoordinatorTools.hpp"

// Project includes
//
#include "QuICC/Timers/StageTimer.hpp"

namespace QuICC {

namespace Transform {

   TransformCoordinatorTools::TransformCoordinatorTools()
   {
   }

   TransformCoordinatorTools::~TransformCoordinatorTools()
   {
   }

   void TransformCoordinatorTools::init(TransformCoordinatorType& rCoord, SharedIForwardGrouper spFwdGrouper, SharedIBackwardGrouper spBwdGrouper, const std::vector<Transform::TransformTree>& forwardTree, const std::vector<Transform::TransformTree>& backwardTree, SharedResolution spRes, const std::map<std::size_t,NonDimensional::SharedINumber>& runOptions)
   {
      StageTimer stage;
      stage.start("initializing transforms");

      // Initialise the transform coordinator
      rCoord.defineTransforms(forwardTree, backwardTree, spRes->sim().spSpatialScheme());
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

      std::vector<ArrayI>  packs;
      // Get the buffer pack sizes for first dimension
      packs.push_back(spFwdGrouper->packs1D(forwardTree));
      packs.push_back(spBwdGrouper->packs1D(backwardTree));

      if(spRes->sim().ss().dimension() == 3)
      {
         // Get the buffer pack sizes for second dimension
         packs.push_back(spFwdGrouper->packs2D(forwardTree));
         packs.push_back(spBwdGrouper->packs2D(backwardTree));
      }

      // Initialise the converters
      rCoord.communicator().initConverter(spRes, packs, spFwdGrouper->split);

      stage.done();
   }
}
}
