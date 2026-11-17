/**
 * @file IVectorEquation.hpp
 * @brief Base for the implementation of a vector equation
 */

#ifndef QUICC_EQUATIONS_IVECTOREQUATION_HPP
#define QUICC_EQUATIONS_IVECTOREQUATION_HPP

// System includes
//
#include <vector>
#include <memory>


// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base for the implementation of a vector equation
    */
   class IVectorEquation: public IFieldEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          */
         explicit IVectorEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          * @param spOptions  Additional options
          */
         explicit IVectorEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IVectorEquation() = default;

         /**
          * @brief Set the smart pointer to the unknown field
          *
          * This is required because the field are not initialised at creation time
          *
          * \param spUnknown Shared pointer to the unknown of the equation
          */
         virtual void setUnknown(Framework::Selector::VariantSharedVectorVariable spUnknown);

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
          * @brief Transfer solver solution to equation unknown
          *
          * @param compId  Component ID
          * @param storage Solver solution
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         template <typename TData> void storeSolution(FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start);

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
          * @brief Get the shared unknown variable
          */
         Framework::Selector::VariantSharedVectorVariable spUnknown() const;

         /**
          * @brief Get backward transform paths
          */
         virtual std::vector<Transform::TransformPath> backwardPaths() override;

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
         Framework::Selector::VariantSharedVectorVariable mspUnknown;
   };

   /// Templated typedef for shared IVectorEquation
   typedef std::shared_ptr<IVectorEquation> SharedIVectorEquation;

   template <typename TData> void IVectorEquation::storeSolution(FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start)
   {
      std::visit([&](auto&& p){this->storeSolutionImpl(p->rDom(0).rPerturbation(), compId, storage, matIdx, start);}, this->spUnknown());
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_IVECTOREQUATION_HPP
