/**
 * @file CartesianExactScalarState.hpp
 * @brief Implementation of the equation to generate exact scalar states in a cartesian geometries
 */

#ifndef QUICC_EQUATIONS_CARTESIANEXACTSCALARSTATE_HPP
#define QUICC_EQUATIONS_CARTESIANEXACTSCALARSTATE_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/SpectralKernels/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the equation to generate exact scalar state in cartesian geometries
    */
   class CartesianExactScalarState: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model Backend
          */
         CartesianExactScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~CartesianExactScalarState() = default;

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
          * @brief Set the spectral state type id
          *
          * @param compId  Field component to act on
          * @param modes   List of spectral modes with amplitude to create
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

   /// Typedef for a shared CartesianExactScalarState
   typedef std::shared_ptr<CartesianExactScalarState> SharedCartesianExactScalarState;

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_CARTESIANEXACTSCALARSTATE_HPP
