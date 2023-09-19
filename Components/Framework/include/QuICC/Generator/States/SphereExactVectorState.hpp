/**
 * @file SphereExactVectorState.hpp
 * @brief Implementation of the equation to generate exact vector states in a sphere
 */

#ifndef QUICC_EQUATIONS_SPHEREEXACTVECTORSTATE_HPP
#define QUICC_EQUATIONS_SPHEREEXACTVECTORSTATE_HPP

// Configuration includes
//

// System includes
//
#include <tuple>
#include <memory>

// External includes
//

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
          * @brief Initialize nonlinear interaction kernel
          *
          * @param force   Force initialization
          */
         virtual void initNLKernel(const bool force = false);

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

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the nonliner integration components
          */
         virtual void setNLComponents();

         /**
          * @brief Transform path tag
          */
         std::size_t mPathTag;

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
