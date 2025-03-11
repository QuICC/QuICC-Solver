/** 
 * @file IModelBuilder.hpp
 * @brief Interface for implementation of a physical model
 */

#ifndef QUICC_MODEL_IMODELBUILDER_HPP
#define QUICC_MODEL_IMODELBUILDER_HPP

// System includes
//
#include <vector>
#include <memory>

// Project includes
//
#include "QuICC/Model/IPhysicalModel.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Interface for the implementation of a physical models
    */
   template <typename TSim> class IModelBuilder: public IPhysicalModel
   {
      public:
         /**
          * @brief Constructor
          */
         IModelBuilder() = default;

         /**
          * @brief Destructor
          */
         virtual ~IModelBuilder() = default;

         /**
          * @brief Add the required equations
          *
          * @param spSim   Shared simulation object
          */
         virtual void addEquations(std::shared_ptr<TSim> spSim) = 0;

         /**
          * @brief Add the required ASCII output files
          *
          * @param spSim   Shared simulation object
          */
         virtual void addAsciiOutputFiles(std::shared_ptr<TSim> spSim) = 0;

         /**
          * @brief Add the required HDF5 output files
          *
          * @param spSim   Shared simulation object
          */
         virtual void addHdf5OutputFiles(std::shared_ptr<TSim> spSim);

         /** 
          * @brief Add the required statistics output files
          * 
          * @param spSim   Shared simulation object
          */
         virtual void addStatsOutputFiles(std::shared_ptr<TSim> spSim);

         /**
          * @brief Set the initial state
          *
          * @param spSim   Shared simulation object
          */
         virtual void setInitialState(std::shared_ptr<TSim> spSim);

      protected:

      private:
   };

   template <typename TSim> void IModelBuilder<TSim>::addHdf5OutputFiles(std::shared_ptr<TSim> spSim)
   {
      // Field IDs iterator
      std::vector<std::size_t> ids = this->backend().fieldIds();

      // Create and add state file to IO
      auto spState = std::make_shared<Io::Variable::StateFileWriter>(spSim->ss().tag(), spSim->ss().has(SpatialScheme::Feature::RegularSpectrum));
      for(auto it = ids.cbegin(); it != ids.cend(); ++it)
      {
         spState->expect(*it);
      }

      // Add extra field names
      ids.clear();
      ids = this->extraFieldIds();
      for(auto it = ids.cbegin(); it != ids.cend(); ++it)
      {
         spState->expect(*it);
      }

      spSim->addHdf5OutputFile(spState);
   }

   template <typename TSim> void IModelBuilder<TSim>::addStatsOutputFiles(std::shared_ptr<TSim>)
   {
   }

   template <typename TSim> void IModelBuilder<TSim>::setInitialState(std::shared_ptr<TSim> spSim)
   {
      // Field IDs iterator
      std::vector<std::size_t> ids = this->backend().fieldIds();

      // Create and add initial state file to IO
      auto spInit = std::make_shared<Io::Variable::StateFileReader>("_initial", spSim->ss().tag(), spSim->ss().has(SpatialScheme::Feature::RegularSpectrum));

      // Set expected field names
      for(auto it = ids.cbegin(); it != ids.cend(); ++it)
      {
         spInit->expect(*it);
      }

      // Add extra field names
      ids.clear();
      ids = this->extraFieldIds();
      for(auto it = ids.cbegin(); it != ids.cend(); ++it)
      {
         spInit->expect(*it);
      }

      // Set simulation state
      spSim->setInitialState(spInit);
   }

}
}

#endif // QUICC_MODEL_IMODELBUILDER_HPP
