/**
 * @file SphereExactVectorState.hpp
 * @brief Implementation of the equation to generate exact vector states in a sphere
 */

#ifndef QUICC_EQUATIONS_SPHEREEXACTVECTORSTATE_HPP
#define QUICC_EQUATIONS_SPHEREEXACTVECTORSTATE_HPP

// System includes
//
#include <tuple>
#include <memory>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/SpectralKernels/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact vector state in a sphere
    */
   class SphereExactVectorState: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model Backend
          */
         SphereExactVectorState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~SphereExactVectorState() = default;

         /**
          * @brief Set backward path ID
          *
          * @param pathId Path ID
          */
         void setBackwardPath(const std::size_t id);

         /**
          * @brief Set forward path ID
          *
          * @param pathId Path ID
          */
         void setForwardPath(const std::size_t id);

         /**
          * @brief Initialize nonlinear interaction kernel
          *
          * @param force   Force initialization
          */
         virtual void initNLKernel(const bool force = false) override;

         /**
          * @brief Set the unknown name and requirements
          *
          * @param name Name of the main output field
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Use nonlinear transform path
          *
          * @param tag ID of the transform path
          */
         void useNonlinearPath(const std::size_t tag);

         /**
          * @brief Set the physical space kernel
          *
          * @param spKernel physical space kernel
          */
         void setPhysicalKernel(const Physical::Kernel::SharedIPhysicalKernel spKernel);

         /**
          * @brief Use noise as physical state
          *
          * @param level Noise level
          */
         void setPhysicalNoise(const MHDFloat level);

         /**
          * @brief Use constant as physical state
          *
          * @param value   Constant value
          */
         void setPhysicalConstant(const MHDFloat value);

         /**
          * @brief Set options for the harmonics states
          *
          * @param compId  ID of the field component
          * @param modes   List of harmonics with amplitude to create
          */
         void setSpectralModes(const FieldComponents::Spectral::Id compId, const Spectral::Kernel::Complex3DMapType& modes);

         /**
          * @brief Initialize constraint kernels
          *
          * @param spMesh  Physical mesh
          */
         virtual void initConstraintKernel(const std::shared_ptr<std::vector<Array> > spMesh) override;

         /**
          * @brief Get backward transform paths
          */
         virtual std::vector<Transform::TransformPath> backwardPaths() override;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements() override;

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling() override;

         /**
          * @brief Set the nonliner integration components
          */
         virtual void setNLComponents() override;

         /**
          * @brief Backward Transform path ID
          */
         std::size_t mBwdPathId;

         /**
          * @brief Forward Transform path ID
          */
         std::size_t mFwdPathId;

      private:
         /**
          * @brief Physical kernel
          */
         Physical::Kernel::SharedIPhysicalKernel mspPhysKernel;
   };

   /// Typedef for a shared SphereExactVectorState
   typedef std::shared_ptr<SphereExactVectorState> SharedSphereExactVectorState;

}
}

#endif // QUICC_EQUATIONS_SPHEREEXACTVECTORSTATE_HPP
