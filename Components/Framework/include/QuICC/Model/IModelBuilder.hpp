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
#include "QuICC/PhysicalNames/Magnetic.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/VelocityZ.hpp"
#include "QuICC/PhysicalNames/Streamfunction.hpp"
#include "QuICC/Diagnostics/InertialWaveCfl.hpp"
#include "QuICC/Diagnostics/TorsionalOscillationCfl.hpp"
#include "QuICC/Diagnostics/ISphericalHydroCfl.hpp"
#include "QuICC/Diagnostics/ISphericalMagneticCfl.hpp"
#include "QuICC/Diagnostics/ShellCfl.hpp"
#include "QuICC/Diagnostics/SphereCfl.hpp"
#include "QuICC/Diagnostics/CartesianCfl.hpp"

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

         /**
          * @brief Add diagnostics
          */
         void addDiagnostics(std::shared_ptr<TSim> spSim);

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

   template <typename TSim> void IModelBuilder<TSim>::addDiagnostics(std::shared_ptr<TSim> spSim)
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

}
}

#endif // QUICC_MODEL_IMODELBUILDER_HPP
