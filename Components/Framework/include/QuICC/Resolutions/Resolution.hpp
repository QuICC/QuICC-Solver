/**
 * @file Resolution.hpp
 * @brief Definition of a resolution object
 */

#ifndef QUICC_RESOLUTION_HPP
#define QUICC_RESOLUTION_HPP

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
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Resolutions/SimulationResolution.hpp"
#include "QuICC/Resolutions/CoreResolution.hpp"
#include "QuICC/Transform/TransformSetup.hpp"

// Forward declarations
namespace QuICC {
   class IndexCounter;

   namespace SpatialScheme {
      class ISpatialScheme;
   }

   namespace Datatypes {
      class ScalarFieldSetup;
   }
}

namespace QuICC {

   /**
    * @brief Definition of a resolution object
    */
   class Resolution
   {
      public:
         /**
          * @brief Constructor
          *
          * @param coreRes    Resolution object for the different CPUs
          * @param simDim     Simulation dimensions
          * @param transDim   Transform dimensions
          */
         Resolution(const std::vector<SharedCoreResolution>& coreRes, const ArrayI& simDim, const ArrayI& transDim);

         /**
          * @brief Empty Destructor
          */
         ~Resolution();

         /**
          * @brief Get the simulation resolution
          */
         SharedCSimulationResolution spSim() const;

         /**
          * @brief Get the simulation resolution
          */
         const SimulationResolution& sim() const;

         /**
          * @brief Get resolution corresponding to local CPU
          */
         SharedCCoreResolution spCpu() const;

         /**
          * @brief Get resolution corresponding to local CPU
          */
         SharedCCoreResolution cpu() const;

         /**
          * @brief Get resolution corresponding to id
          *
          * @param id CPU id
          */
         SharedCCoreResolution cpu(const int id) const;

         /**
          * @brief Get Number of CPUs
          */
         int nCpu() const;

         /**
          * @brief Get the transform setup
          *
          * @param dim Dimension for which to get setup
          */
         Transform::SharedTransformSetup spTransformSetup(const Dimensions::Transform::Id dim) const;

         /**
          * @brief Add a transform setup
          *
          * @param dim     Dimension corresponding to setup
          * @param spSetup Transform setup
          */
         void addTransformSetup(const Dimensions::Transform::Id dim, Transform::SharedTransformSetup spSetup);

         /**
          * @brief Get the forward scalar field setup
          *
          * @param dim Dimension for which to get setup
          */
         std::shared_ptr<Datatypes::ScalarFieldSetup> spFwdSetup(const Dimensions::Transform::Id dim) const;

         /**
          * @brief Get the backward scalar field setup
          *
          * @param dim Dimension for which to get setup
          */
         std::shared_ptr<Datatypes::ScalarFieldSetup> spBwdSetup(const Dimensions::Transform::Id dim) const;

         /**
          * @brief Get the forward scalar field setup for last dimension
          */
         std::shared_ptr<Datatypes::ScalarFieldSetup> spPhysicalSetup() const;

         /**
          * @brief Get the backward scalar field setup for the first dimension
          */
         std::shared_ptr<Datatypes::ScalarFieldSetup> spSpectralSetup() const;

         /**
          * @brief Set the box scale for the periodic box dimensions
          */
         void setBoxScale(const Array& boxScale);

         /**
          * @brief Set the spatial scheme
          */
         void setSpatialScheme(std::shared_ptr<SpatialScheme::ISpatialScheme> spScheme);

         /**
          * @brief Set the index counter
          */
         void setIndexCounter(std::shared_ptr<IndexCounter> spCounter);

         /**
          * @brief Get the index counter
          */
         const IndexCounter& counter() const;

         /**
          * @brief Build restriction for MPI sparse solver
          */
         void buildRestriction(std::vector<int>& rSlow, std::vector<std::vector<int> >& rMiddle, const int k);

      protected:

      private:
         /**
          * @brief Extract local ID from core resolution
          */
         void setLocalId();

         /**
          * @brief Initialise the simulation resolution
          *
          * @param simDim     Simulation dimensions
          * @param transDim   Transform dimensions
          */
         void initSimResolution(const ArrayI& simDim, const ArrayI& transDim);

         /**
          * @brief Local ID
          */
         int mLocalId;

         /**
          * @brief Shared simulation resolution
          */
         SharedSimulationResolution mspSim;

         /**
          * @brief Storage for all the core resolutions
          */
         std::vector<SharedCoreResolution>   mCores;

         /**
          * @brief Storage for the transform setups
          */
         std::map<Dimensions::Transform::Id,Transform::SharedTransformSetup>   mTSetups;

         /**
          * @brief Storage for the index counter
          */
         std::shared_ptr<IndexCounter>   mspCounter;
   };

   /// Typedef for a shared pointer to a Resolution object
   typedef std::shared_ptr<Resolution>   SharedResolution;

   /// Typedef for a shared pointer to a const Resolution object
   typedef std::shared_ptr<const Resolution>   SharedCResolution;

}

#endif // QUICC_RESOLUTION_HPP
