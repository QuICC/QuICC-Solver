/**
 * @file ShellExactScalarState.hpp
 * @brief Implementation of the equation to generate exact scalar states in a spherical shell
 */

#ifndef QUICC_EQUATIONS_SHELLEXACTSCALARSTATE_HPP
#define QUICC_EQUATIONS_SHELLEXACTSCALARSTATE_HPP

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
#include "QuICC/Typedefs.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/SpectralKernels/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact scalar state in a spherical shell
    */
   class ShellExactScalarState: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         ShellExactScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ShellExactScalarState();

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

   /// Typedef for a shared ShellExactScalarState
   typedef std::shared_ptr<ShellExactScalarState> SharedShellExactScalarState;

}
}

#endif // QUICC_EQUATIONS_SHELLEXACTSCALARSTATE_HPP
