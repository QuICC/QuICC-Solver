/** 
 * @file IPhysicalModel.hpp
 * @brief Interface for implementation of a physical model
 */

#ifndef QUICC_MODEL_IPHYSICALMODEL_HPP
#define QUICC_MODEL_IPHYSICALMODEL_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <vector>
#include <set>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Enums/VectorFormulation.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Arithmetics/registerAll.hpp"
#include "QuICC/ModelOperator/registerAll.hpp"
#include "QuICC/ModelOperatorBoundary/registerAll.hpp"
#include "QuICC/NonDimensional/registerAll.hpp"
#include "QuICC/PhysicalNames/registerAll.hpp"
#include "QuICC/RuntimeStatus/registerAll.hpp"
#include "QuICC/SolveTiming/registerAll.hpp"
#include "QuICC/PseudospectralTag/registerAll.hpp"

namespace QuICC {

namespace Model {

   /**
    * @brief Interface for the implementation of a physical models
    */
   template <typename TSim, typename TState, typename TVis> class IPhysicalModel
   {
      public:
         /**
          * @brief Constructor
          */
         IPhysicalModel() = default;

         /**
          * @brief Destructor
          */
         virtual ~IPhysicalModel() = default;

         /**
          * @brief Initialize model
          */
         virtual void init();

         /// Formulation used for vector fields
         virtual VectorFormulation::Id SchemeFormulation() = 0;

         /**
          * @brief Version string of model
          */
         virtual std::string version() const = 0;

         /**
          * @brief XML configuration tags for model
          */
         virtual std::map<std::string,std::map<std::string,int> > configTags() const;

         /**
          * @brief Configure additional features set at run time
          */
         virtual void configure(const std::set<SpatialScheme::Feature>& f);

         /**
          * @brief Add extra field IDs (example: imposed fields)
          */
         virtual std::vector<std::size_t> extraFieldIds() const;

         /**
          * @brief Add the required equations
          *
          * @param spSim   Shared simulation object
          */
         virtual void addEquations(std::shared_ptr<TSim> spSim) = 0;

         /**
          * @brief Add the initial state generation equations
          *
          * @param spGen   Shared generator object
          */
         virtual void addStates(std::shared_ptr<TState> spGen) = 0;

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

         /**
          * @brief Get model generator
          */
         const IModelBackend& backend() const;

         /**
          * @brief Get model generator
          */
         std::shared_ptr<IModelBackend> spBackend() const;

         /**
          * @brief Interface to adding ASCII output file
          */
         template <typename T> void enableAsciiFile(const std::string tag, const std::string prefix, const std::size_t id, std::shared_ptr<TSim> spSim);

      protected:
         /**
          * @brief Register Named IDs needed for simulation
          */
         virtual void registerNames();

         /**
          * @brief Model generator
          */
         std::shared_ptr<IModelBackend> mpBackend;

      private:
   };

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::init()
   {
      this->registerNames();
   }

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::registerNames()
   {
      // Arithmetics names
      Arithmetics::registerAll();
      // ModelOperator names
      ModelOperator::registerAll();
      // ModelOperatorBoundary names
      ModelOperatorBoundary::registerAll();
      // NonDimensional names
      NonDimensional::registerAll();
      // Physical names
      PhysicalNames::registerAll();
      // RuntimeStatus names
      RuntimeStatus::registerAll();
      // SolveTiming names
      SolveTiming::registerAll();
      // PseudospectralTag names
      PseudospectralTag::registerAll();
   }

   template <typename TSim, typename TState, typename TVis> std::vector<std::size_t> IPhysicalModel<TSim,TState,TVis>::extraFieldIds() const
   {
      std::vector<std::size_t> extra;

      return extra;
   }

   template <typename TSim, typename TState, typename TVis> std::map<std::string, std::map<std::string,int> > IPhysicalModel<TSim,TState,TVis>::configTags() const
   {
      std::map<std::string, std::map<std::string,int> > tags;

      return tags;
   }

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::configure(const std::set<SpatialScheme::Feature>& f)
   {
      this->mpBackend->enableGalerkin(f.count(SpatialScheme::Feature::GalerkinBasis));
   }

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::setVisualizationState(std::shared_ptr<TVis> spVis)
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

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::addHdf5OutputFiles(std::shared_ptr<TSim> spSim)
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

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::addStatsOutputFiles(std::shared_ptr<TSim>)
   {
   }

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::setInitialState(std::shared_ptr<TSim> spSim)
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

   template <typename TSim, typename TState, typename TVis> const IModelBackend& IPhysicalModel<TSim,TState,TVis>::backend() const
   {
      return *this->mpBackend;
   }

   template <typename TSim, typename TState, typename TVis> std::shared_ptr<IModelBackend> IPhysicalModel<TSim,TState,TVis>::spBackend() const
   {
      return this->mpBackend;
   }

   template <typename TSim, typename TState, typename TVis> template <typename T> void IPhysicalModel<TSim,TState,TVis>::enableAsciiFile(const std::string tag, const std::string prefix, const std::size_t id, std::shared_ptr<TSim> spSim)
   {
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<T>(prefix, spSim->ss().tag());
         spFile->expect(id);
         if((spSim->config().model(tag).count("numbered") > 0) && spSim->config().model(tag).at("numbered"))
         {
            spFile->numberOutput();
         }
         if(spSim->config().model(tag).count("only_every") > 0)
         {
            spFile->onlyEvery(spSim->config().model(tag).at("only_every"));
         }
         spSim->addAsciiOutputFile(spFile);
      }
   }

}
}

#endif // QUICC_MODEL_IPHYSICALMODEL_HPP
