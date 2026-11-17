/**
 *
 * @file IScalarEquation.hpp
 * @brief Base for the implementation of a scalar equation
 */

#ifndef QUICC_EQUATIONS_ISCALAREQUATION_HPP
#define QUICC_EQUATIONS_ISCALAREQUATION_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/CouplingFeature.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base for the implementation of a scalar equation
    */
   class IScalarEquation: public IFieldEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          */
         explicit IScalarEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          * @param spOptions  Additional options
          */
         explicit IScalarEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IScalarEquation() = default;

         /**
          * @brief Set the shared pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         virtual void setUnknown(Framework::Selector::VariantSharedScalarVariable spUnknown);

         /**
          * @brief Access the shared resolution
          */
         virtual SharedResolution spRes() const override;

         /**
          * @brief Access the resolution
          */
         virtual const Resolution& res() const override;

         /**
          * @brief Get the number of spectral components
          */
         int nSpectral() const final;

         /**
          * @brief Get vector spectral component range
          */
         SpectralComponent_range spectralRange() const final;

         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices() override;

         /**
          * @brief Generic model operator dispatcher to python scripts
          */
         virtual void buildModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id comp, const int matIdx, const std::size_t bcType) const override;

         /**
          * @brief Get the shared pointer to unknown variable
          */
         Framework::Selector::VariantSharedScalarVariable spUnknown() const;

         /**
          * @brief Get backward transform paths
          */
         virtual std::vector<Transform::TransformPath> backwardPaths() override;

         /**
          * @brief Set spectral constraint kernel overload for scalar case
          */
         void setConstraintKernel(Spectral::Kernel::SharedISpectralKernel spKernel);
         using IFieldEquation::setConstraintKernel;

         /**
          * @brief Set spectral source kernel overload for scalar case
          */
         virtual void setSrcKernel(Spectral::Kernel::SharedISpectralKernel spKernel);
         using IFieldEquation::setSrcKernel;

         /**
          * @brief Set unknown field to bad value
          *
          * @param compId  Component ID
          */
         void corruptUnknown(FieldComponents::Spectral::Id compId);

      protected:
         /**
          * @brief Get backward transform paths
          *
          * @param pathId  ID of enabled path
          */
         std::vector<Transform::TransformPath> defaultBackwardPaths(const std::size_t pathId) const;

         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents() override;

         /**
          * @brief Get backward transform paths
          */
         virtual std::vector<bool> disabledBackwardPaths() const override;

         /**
          * @brief Set the galerkin stencil
          */
         virtual void setGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const override;

         /**
          * @brief Set the explicit matrix operator
          */
         virtual void setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId, const SpectralFieldId fieldId, const int matIdx) const override;

         /**
          * @brief Build coupling information from Python scripts
          */
         void defineCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const std::map<CouplingFeature,bool>& features);

      private:
         /**
          * @brief Storage for the shared scalar variable
          */
         Framework::Selector::VariantSharedScalarVariable mspUnknown;
   };

   /// Typedef for shared IScalarEquation
   typedef std::shared_ptr<IScalarEquation> SharedIScalarEquation;

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_ISCALAREQUATION_HPP
