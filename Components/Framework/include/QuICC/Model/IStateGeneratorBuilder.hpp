/** 
 * @file IStateGeneratorBuilder.hpp
 * @brief Interface for a state generator builder
 */

#ifndef QUICC_MODEL_ISTATEGENERATORBUILDER_HPP
#define QUICC_MODEL_ISTATEGENERATORBUILDER_HPP

// System includes
//
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
   template <typename TState> class IStateGeneratorBuilder: public IPhysicalModel
   {
      public:
         /**
          * @brief Constructor
          */
         IStateGeneratorBuilder() = default;

         /**
          * @brief Destructor
          */
         virtual ~IStateGeneratorBuilder() = default;

         /**
          * @brief Add the initial state generation equations
          *
          * @param spGen   Shared generator object
          */
         virtual void addStates(std::shared_ptr<TState> spGen) = 0;

         /**
          * @brief Set the state generator initial state
          *
          * @param spGen   Shared generator
          */
         void setGeneratorState(std::shared_ptr<TState> spGen);

      protected:
         /**
          * @brief Set default state generator initial state
          *
          * @param spGen   Shared generator
          */
         void setDefaultGeneratorState(std::shared_ptr<TState> spGen);

      private:
   };

   template <typename TState> void IStateGeneratorBuilder<TState>::setGeneratorState(std::shared_ptr<TState> spGen)
   {
   }

   template <typename TState> void IStateGeneratorBuilder<TState>::setDefaultGeneratorState(std::shared_ptr<TState> spGen)
   {
      // Field IDs iterator
      std::vector<std::size_t> ids = this->backend().fieldIds();

      // Create and add initial state file to IO
      auto spIn = std::make_shared<Io::Variable::StateFileReader>("_initial", spGen->ss().tag(), spGen->ss().has(SpatialScheme::Feature::RegularSpectrum));

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
      spGen->setInitialState(spIn);
   }

}
}

#endif // QUICC_MODEL_ISTATEGENERATORBUILDER_HPP
