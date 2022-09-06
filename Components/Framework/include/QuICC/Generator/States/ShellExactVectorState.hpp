/**
 * @file ShellExactVectorState.hpp
 * @brief Implementation of the equation to generate exact vector states in a shell
 */

#ifndef QUICC_EQUATIONS_SHELLEXACTVECTORSTATE_HPP
#define QUICC_EQUATIONS_SHELLEXACTVECTORSTATE_HPP

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
#include "QuICC/Equations/IVectorEquation.hpp"
#include "QuICC/SpectralKernels/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact vector state in a shell
    */
   class ShellExactVectorState: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         ShellExactVectorState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ShellExactVectorState();

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
          * @brief Set options for the harmonics states
          *
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

      private:
         /**
          * @brief Physical kernel
          */
         Physical::Kernel::SharedIPhysicalKernel mspPhysKernel;
   };

   /// Typedef for a shared ShellExactVectorState
   typedef std::shared_ptr<ShellExactVectorState> SharedShellExactVectorState;

}
}

#endif // QUICC_EQUATIONS_SHELLEXACTVECTORSTATE_HPP
