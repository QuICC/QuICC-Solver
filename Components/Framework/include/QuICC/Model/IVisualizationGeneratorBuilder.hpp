/** 
 * @file IVisualizationGeneratorBuilder.hpp
 * @brief Interface for implementation of a physical model
 */

#ifndef QUICC_MODEL_IVISUALIZATIONGENERATORBUILDER_HPP
#define QUICC_MODEL_IVISUALIZATIONGENERATORBUILDER_HPP

// System includes
//
#include <vector>
#include <memory>

// Project includes
//
#include "QuICC/Model/IPhysicalModel.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Interface for the implementation of a physical models
    */
   template <typename TVis> class IVisualizationGeneratorBuilder: public IPhysicalModel
   {
      public:
         /**
          * @brief Constructor
          */
         IVisualizationGeneratorBuilder() = default;

         /**
          * @brief Destructor
          */
         virtual ~IVisualizationGeneratorBuilder() = default;

         /**
          * @brief Add the visualization generation equations
          *
          * @param spGen   Shared visualization generator
          */
         virtual void addVisualizers(std::shared_ptr<TVis> spVis) = 0;

         /**
          * @brief Set the visualization initial state
          *
          * @param spSim   Shared visualization generator
          */
         virtual void setVisualizationState(std::shared_ptr<TVis> spVis);

      protected:

      private:
   };

   template <typename TVis> void IVisualizationGeneratorBuilder<TVis>::setVisualizationState(std::shared_ptr<TVis> spVis)
   {
      // Field IDs iterator
      std::vector<std::size_t> ids = this->backend().fieldIds();

      // Create and add initial state file to IO
      auto spIn = std::make_shared<Io::Variable::StateFileReader>("4Visu", spVis->ss().tag(), spVis->ss().has(SpatialScheme::Feature::RegularSpectrum));

      // Set expected field names
      for(auto it = ids.cbegin(); it != ids.cend(); ++it)
      {
         spIn->expect(*it);
      }

      // Add extra field names
      ids.clear();
      ids = this->extraFieldIds();
      for(auto it = ids.cbegin(); it != ids.cend(); ++it)
      {
         spIn->expect(*it);
      }

      // Set simulation state
      spVis->setInitialState(spIn);
   }

}
}

#endif // QUICC_MODEL_IVISUALIZATIONGENERATORBUILDER_HPP
