/**
 *
 * @file IScalarEquation.hpp
 * @brief Base for the implementation of a scalar equation
 */

#ifndef QUICC_EQUATIONS_ISCALAREQUATION_HPP
#define QUICC_EQUATIONS_ISCALAREQUATION_HPP

// Configuration includes
//

// System includes
//
#include <memory>
#include <stdexcept>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/DecoupledComplexInternal.hpp"
#include "QuICC/Equations/CouplingFeature.hpp"
#include "QuICC/Solver/SparseSolver.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolverTools.hpp"

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
          * \param spEqParams Shared equation parameters
          */
         explicit IScalarEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IScalarEquation();

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

      protected:
         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents() override;

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

   /**
    * @brief Solve for galerkin unknown using the stencil
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    */
   template <typename TData> void solveStencilUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Transfer nonlinear spectral values from unknown to solver
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData> void copyNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   template <typename TData> void IScalarEquation::storeSolution(FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);

      std::visit([&](auto&& p){this->storeSolutionImpl(p->rDom(0).rPerturbation(), compId, storage, matIdx, start);}, this->spUnknown());
   }

   template <typename TData> void solveStencilUnknown(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Create temporary storage for tau data
      TData tmp(eq.couplingInfo(compId).tauN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      std::visit(
            [&](auto&& p)
            {
            Equations::copyUnknown(eq, p->dom(0).perturbation(), compId, tmp, matIdx, 0, false, true);
            }, eq.spUnknown());
      TData rhs(eq.couplingInfo(compId).galerkinN(matIdx), eq.couplingInfo(compId).rhsCols(matIdx));
      if(eq.res().sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
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

      // Check that square stencil was generated (Python setup is wrong if matrix is not square)
      if(stencil.rows() != stencil.cols())
      {
      	throw std::logic_error("Stencil setup is wrong and did not produce a square matrix");
      }

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

   template <typename TData> void copyNonlinear(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);
      assert((!eq.couplingInfo(compId).isGalerkin() || eq.couplingInfo(compId).indexType() != CouplingIndexType::SINGLE) && "Current version does not support galerkin basis");

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

      /// Nonlinear computation took place but no quasi-inverse is required
      } else if(eq.couplingInfo(compId).hasNonlinear())
      {
         // simply copy values from unknown
         std::visit([&](auto&& p){copyUnknown(eq, p->dom(0).perturbation(), compId, storage, matIdx, start, true, false);}, eq.spUnknown());
      }
   }
}
}

#endif // QUICC_EQUATIONS_ISCALAREQUATION_HPP
