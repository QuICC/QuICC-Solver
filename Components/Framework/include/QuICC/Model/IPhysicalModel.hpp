/** 
 * @file IPhysicalModel.hpp
 * @brief Interface for implementation of a physical model
 */

#ifndef QUICC_MODEL_IPHYSICALMODEL_HPP
#define QUICC_MODEL_IPHYSICALMODEL_HPP

// System includes
//
#include <string>
#include <vector>
#include <set>
#include <memory>

// Project includes
//
#include "QuICC/Enums/VectorFormulation.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"
#include "QuICC/Arithmetics/registerAll.hpp"
#include "QuICC/ModelOperator/registerAll.hpp"
#include "QuICC/ModelOperatorBoundary/registerAll.hpp"
#include "QuICC/NonDimensional/registerAll.hpp"
#include "QuICC/PhysicalNames/Streamfunction.hpp"
#include "QuICC/PhysicalNames/registerAll.hpp"
#include "QuICC/RuntimeStatus/registerAll.hpp"
#include "QuICC/PseudospectralTag/registerAll.hpp"
#include "QuICC/Transform/Reductor/registerAll.hpp"
#include "QuICC/SolveTiming/registerAll.hpp"
#include "QuICC/Diagnostics/CartesianCfl.hpp"
#include "QuICC/Diagnostics/ShellCfl.hpp"
#include "QuICC/Diagnostics/SphereCfl.hpp"
#include "QuICC/Diagnostics/InertialWaveCfl.hpp"
#include "QuICC/Diagnostics/TorsionalOscillationCfl.hpp"
#include "QuICC/Diagnostics/ISphericalHydroCfl.hpp"
#include "QuICC/Diagnostics/ISphericalMagneticCfl.hpp"

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
          * @brief Tune the spatial scheme (for example change mesher)
          *
          * @param spScheme   Spatial scheme
          */
         virtual void tuneScheme(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme) {};

         /**
          * @brief Initialize model
          */
         virtual void init();

         /**
          * @brief Formulation used for vector fields
          */
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
          * @brief Set the state generator initial state
          *
          * @param spGen   Shared generator object
          */
         virtual void setGeneratorState(std::shared_ptr<TState> spGen);

         /**
          * @brief Add the visualization generation equations
          *
          * @param spGen   Shared visualization generator
          */
         virtual void addVisualizers(std::shared_ptr<TVis> spVis) = 0;

         /**
          * @brief Set the visualization initial state
          *
          * @param spVis   Shared visualization generator
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
          * @brief Add diagnostic
          *
          * Default implementation provides standard CFL calculation
          *
          * @param spSim   Shared simulation object
          */
         virtual void addDiagnostics(std::shared_ptr<TSim> spSim);

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
         template <typename T> std::shared_ptr<T> enableAsciiFile(const std::string tag, const std::string prefix, const std::size_t id, std::shared_ptr<TSim> spSim);

      protected:
         /**
          * @brief Set default generator initial state
          *
          * @param spGen   Shared generator object
          */
         void setDefaultGeneratorState(std::shared_ptr<TState> spGen);

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
      Transform::Reductor::registerAll();
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
      // Propagate Galerkin flag
      this->mpBackend->enableGalerkin(f.count(SpatialScheme::Feature::GalerkinBasis));

      // Propagate split 4th order equations flag
      this->mpBackend->enableSplitEquation(f.count(SpatialScheme::Feature::SplitFourthOrder));
   }

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::addDiagnostics(std::shared_ptr<TSim> spSim)
   {
      const auto& params = spSim->eqParams()->map();

      auto fields = this->backend().fieldIds();

      // Create a toroidal/poloidal spherical shell default CFL
      if (spSim->ss().has(SpatialScheme::Feature::ShellGeometry) &&
          spSim->ss().formulation() == VectorFormulation::TORPOL)
      {
         const MHDFloat courant = 0.65;

         auto iwCfl = std::make_shared<Diagnostics::InertialWaveCfl>(params, courant);
         spSim->addCfl(iwCfl);

         auto toCfl = std::make_shared<Diagnostics::TorsionalOscillationCfl>(params, courant);
         spSim->addCfl(toCfl);

         bool hasVel = std::find(fields.begin(), fields.end(), PhysicalNames::Velocity::id()) != fields.end();
         bool hasMag = std::find(fields.begin(), fields.end(), PhysicalNames::Magnetic::id()) != fields.end();

         if(hasMag && hasVel)
         {
            auto locCfl = std::make_shared<Diagnostics::ShellCfl<Diagnostics::ISphericalMagneticCfl>>(params, courant);
            locCfl->defineVelocity(PhysicalNames::Velocity::id());
            locCfl->defineMagnetic(PhysicalNames::Magnetic::id());
            spSim->addCfl(locCfl);
         }
         else if(hasVel)
         {
            auto locCfl = std::make_shared<Diagnostics::ShellCfl<Diagnostics::ISphericalHydroCfl>>(params, courant);
            locCfl->defineVelocity(PhysicalNames::Velocity::id());
            spSim->addCfl(locCfl);
         }
      }
      // Create a toroidal/poloidal full sphere default CFL
      else if (spSim->ss().has(SpatialScheme::Feature::SphereGeometry) &&
               spSim->ss().formulation() == VectorFormulation::TORPOL)
      {
         const MHDFloat courant = 0.65;

         auto iwCfl = std::make_shared<Diagnostics::InertialWaveCfl>(params, courant);
         spSim->addCfl(iwCfl);

         auto toCfl = std::make_shared<Diagnostics::TorsionalOscillationCfl>(params, courant);
         spSim->addCfl(toCfl);

         bool hasVel = std::find(fields.begin(), fields.end(), PhysicalNames::Velocity::id()) != fields.end();
         bool hasMag = std::find(fields.begin(), fields.end(), PhysicalNames::Magnetic::id()) != fields.end();

         if(hasMag && hasVel)
         {
            auto locCfl = std::make_shared<Diagnostics::SphereCfl<Diagnostics::ISphericalMagneticCfl>>(params, courant);
            locCfl->defineVelocity(PhysicalNames::Velocity::id());
            locCfl->defineMagnetic(PhysicalNames::Magnetic::id());
            spSim->addCfl(locCfl);
         }
         else if(hasVel)
         {
            auto locCfl = std::make_shared<Diagnostics::SphereCfl<Diagnostics::ISphericalHydroCfl>>(params, courant);
            locCfl->defineVelocity(PhysicalNames::Velocity::id());
            spSim->addCfl(locCfl);
         }
      }
      // Create a toroidal/poloidal cartesian default CFL
      else if (spSim->ss().has(SpatialScheme::Feature::CartesianGeometry) &&
          spSim->ss().formulation() == VectorFormulation::TORPOL)
      {
         const MHDFloat courant = 0.65;

         bool hasVel = std::find(fields.begin(), fields.end(), PhysicalNames::Velocity::id()) != fields.end();

         if(hasVel)
         {
            auto locCfl = std::make_shared<Diagnostics::CartesianCfl>(params, courant);
            locCfl->defineVelocity(PhysicalNames::Velocity::id());
            spSim->addCfl(locCfl);
         }
      }
      if (spSim->ss().has(SpatialScheme::Feature::CartesianGeometry))
      {
         bool hasStream = std::find(fields.begin(), fields.end(), PhysicalNames::Streamfunction::id()) != fields.end();
         bool hasVelZ = std::find(fields.begin(), fields.end(), PhysicalNames::VelocityZ::id()) != fields.end();

         if(hasVelZ && hasStream)
         {
            auto locCfl = std::make_shared<Diagnostics::CartesianCfl>(params, 0.65);
            locCfl->defineVelocity(PhysicalNames::Streamfunction::id());
            spSim->addCfl(locCfl);
         }
      }
   }

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::setGeneratorState(std::shared_ptr<TState> spGen)
   {
   }

   template <typename TSim, typename TState, typename TVis> void IPhysicalModel<TSim,TState,TVis>::setDefaultGeneratorState(std::shared_ptr<TState> spGen)
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

   template <typename TSim, typename TState, typename TVis> template <typename T> std::shared_ptr<T> IPhysicalModel<TSim,TState,TVis>::enableAsciiFile(const std::string tag, const std::string prefix, const std::size_t id, std::shared_ptr<TSim> spSim)
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

         return spFile;
      }
      else
      {
         return nullptr;
      }
   }

}
}

#endif // QUICC_MODEL_IPHYSICALMODEL_HPP
