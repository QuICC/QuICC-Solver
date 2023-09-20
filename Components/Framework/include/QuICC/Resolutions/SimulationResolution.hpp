/**
 * @file SimulationResolution.hpp
 * @brief Definition of a simulation resolution object
 */

#ifndef QUICC_SIMULATIONRESOLUTION_HPP
#define QUICC_SIMULATIONRESOLUTION_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <memory>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"

// Forward declarations
namespace QuICC {
   namespace SpatialScheme {
      class ISpatialScheme;
   }
}

namespace QuICC {

   /**
    * @brief Definition of a simulation resolution object
    */
   class SimulationResolution
   {
      public:
         /**
          * @brief Constructor
          *
          * @param phys    Dimensions for the physical space
          * @param spec    Dimensions for the spectral space
          * @param trans    Dimensions for the transform space
          */
         SimulationResolution(const ArrayI& phys, const ArrayI& spec, const ArrayI& trans);

         /**
          * @brief Empty Destructor
          */
         ~SimulationResolution();

         /**
          * @brief Get simulation's dimensions
          *
          * @param simId  ID of the simulation dimension (SIM1D, SIM2D, SIM3D)
          * @param spaceId ID of the space (PHYSICAL, SPECTRAL)
          */
         int dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Get the array of dimensions 1D, 2D, 3D
          *
          * @param spaceId ID of the space (PHYSICAL, SPECTRAL)
          */
         const ArrayI& dimensions(const Dimensions::Space::Id spaceId) const;

         /**
          * @brief Get the box scale
          *
          * @para id ID of the simulation dimension
          */
         MHDFloat boxScale(const Dimensions::Simulation::Id id) const;

         /**
          * @brief Set the box scale (if not it will be initialised to 1)
          */
         void setBoxScale(const Array& boxScale);

         /**
          * @brief Set the spatial scheme
          */
         void setSpatialScheme(std::shared_ptr<SpatialScheme::ISpatialScheme>spScheme);

         /**
          * @brief Get spatial scheme
          */
         const SpatialScheme::ISpatialScheme& ss() const;

         /**
          * @brief Get spatial scheme
          */
         std::shared_ptr<const SpatialScheme::ISpatialScheme> spSpatialScheme() const;

      protected:

      private:
         /**
          * brief Storage
          */
         std::map<Dimensions::Space::Id, ArrayI> mDim;

         /**
          * @brief Storage for the box scale
          */
         Array mBoxScale;

         /**
          * @brief Storage for the spatial scheme
          */
         std::shared_ptr<SpatialScheme::ISpatialScheme>  mspSpatialScheme;
   };


   /// Typedef for a shared pointer to a SimulationResolution object
   typedef std::shared_ptr<SimulationResolution>   SharedSimulationResolution;

   /// Typedef for a shared pointer to a const SimulationResolution object
   typedef std::shared_ptr<const SimulationResolution>   SharedCSimulationResolution;

}

#endif // QUICC_SIMULATIONRESOLUTION_HPP
