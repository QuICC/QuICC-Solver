/**
 * @file IVectorEquation.hpp
 * @brief Base for the implementation of a vector equation
 */

#ifndef QUICC_EQUATIONS_IVECTOREQUATION_HPP
#define QUICC_EQUATIONS_IVECTOREQUATION_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <memory>
#include <stdexcept>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/DecoupledComplexInternal.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolverTools.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base for the implementation of a vector equation
    */
   class IVectorEquation: public IFieldEquation
   {
      public:
         /// Typedef for the the spectral field component ID iterator
         typedef std::vector<FieldComponents::Spectral::Id>::const_iterator   SpectralComponent_iterator;

         /// Typedef for the the spectral field component ID iterator range
         typedef std::pair<SpectralComponent_iterator,SpectralComponent_iterator>  SpectralComponent_range;

         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IVectorEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IVectorEquation();

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
         int nSpectral() const;

         /**
          * @brief Get vector spectral component range
          */
         SpectralComponent_range spectralRange() const;

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

      protected:
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

   /**
    * @brief Solve for galerkin unknown using the stencil
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    */
   template <typename TData> void solveStencilUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Transfer nonlinear spectral values from unknown to solver
    *
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void copyNonlinear(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   template <typename TData> void IVectorEquation::storeSolution(FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start)
   {
      std::visit([&](auto&& p){this->storeSolutionImpl(p->rDom(0).rPerturbation(), compId, storage, matIdx, start);}, this->spUnknown());
   }

   template <typename TData> void solveStencilUnknown(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Create temporary storage for tau data
      TData tmp(eq.couplingInfo(compId).tauN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      std::visit([&](auto&& p){Equations::copyUnknown(eq, p->dom(0).perturbation(), compId, tmp, matIdx, 0, false, true);}, eq.spUnknown());
      TData rhs(eq.couplingInfo(compId).galerkinN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      if(eq.res().sim().ss().has(SpatialScheme::Feature::CylinderGeometry))
      {
         Datatypes::internal::setTopBlock(rhs, 0, eq.couplingInfo(compId).galerkinN(matIdx), eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL), eq.couplingInfo(compId).galerkinShift(matIdx, 0), tmp);
      } else
      {
         Datatypes::internal::setTopBlock(rhs, 0, eq.couplingInfo(compId).galerkinN(matIdx), tmp);
      }

      // Get a restricted stencil matrix
      SparseMatrix stencil(eq.couplingInfo(compId).galerkinN(matIdx),eq.couplingInfo(compId).galerkinN(matIdx));
      eq.dispatchGalerkinStencil(compId, stencil, matIdx, eq.res(), eq.couplingInfo(compId).couplingTools().getIndexes(eq.res(), matIdx), true);
      stencil.makeCompressed();

      // Create solver and factorize stencil
      Framework::Selector::SparseSolver<SparseMatrix> solver;
      solver.compute(stencil);
      // Safety assert for successful factorisation
      if(solver.info() != Eigen::Success)
      {
         throw std::logic_error("Stencil factorization for initial solution failed!");
      }

      // solve for galerkin expansion
      TData lhs(eq.couplingInfo(compId).galerkinN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      Solver::internal::solveWrapper(lhs, solver, rhs);
      Datatypes::internal::setTopBlock(storage, start, eq.couplingInfo(compId).galerkinN(matIdx), lhs);
   }

   template <typename TData> void copyNonlinear(const IVectorEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Check if a nonlinear computation took place and a quasi-inverse has to be applied
      if(eq.couplingInfo(compId).hasNonlinear() && eq.couplingInfo(compId).hasQuasiInverse())
      {
         // Temporary storage is required
         TData tmp;
         tmp = TData(eq.couplingInfo(compId).tauN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));

         // simply copy values from unknown
         std::visit([&](auto&& p){copyUnknown(eq, p->dom(0).perturbation(), compId, tmp, matIdx, 0, false, true);}, eq.spUnknown());

         // Multiply nonlinear term by quasi-inverse
         applyQuasiInverse(eq, compId, storage, start, matIdx, 0, tmp);

      /// Nonlinear computation took place but no quas-inverse is required
      } else if(eq.couplingInfo(compId).hasNonlinear())
      {
         // simply copy values from unknown
         std::visit([&](auto&& p){copyUnknown(eq, p->dom(0).perturbation(), compId, storage, matIdx, start, true, false);}, eq.spUnknown());
      }
   }

}
}

#endif // QUICC_EQUATIONS_IVECTOREQUATION_HPP
