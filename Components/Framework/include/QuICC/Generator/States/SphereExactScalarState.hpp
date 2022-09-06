/**
 * @file SphereExactScalarState.hpp
 * @brief Implementation of the equation to generate exact scalar states in a sphere
 */

#ifndef QUICC_EQUATIONS_SPHEREEXACTSCALARSTATE_HPP
#define QUICC_EQUATIONS_SPHEREEXACTSCALARSTATE_HPP

// Configuration includes
//

// System includes
//
#include <map>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/SpectralKernels/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact scalar state in a sphere
    */
   class SphereExactScalarState: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         SphereExactScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~SphereExactScalarState();

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false);

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const std::size_t name);

         /**
          * @brief Set the physical space kernel
          */
         void setPhysicalKernel(const Physical::Kernel::SharedIPhysicalKernel spKernel);

         /**
          * @brief Use noise as physical state
          */
         void setPhysicalNoise(const MHDFloat level);

         /**
          * @brief Use constant as physical state
          */
         void setPhysicalConstant(const MHDFloat value);

         /**
          * @brief Set the spectral state type id
          */
         void setSpectralModes(const Spectral::Kernel::Complex3DMapType& modes);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
         /**
          * @brief Physical kernel
          */
         Physical::Kernel::SharedIPhysicalKernel mspPhysKernel;
   };

   /// Typedef for a shared SphereExactScalarState
   typedef std::shared_ptr<SphereExactScalarState> SharedSphereExactScalarState;

}
}

#endif // QUICC_EQUATIONS_SPHEREEXACTSCALARSTATE_HPP
