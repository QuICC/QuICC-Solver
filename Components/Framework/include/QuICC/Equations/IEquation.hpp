/**
 * @file IEquation.hpp
 * @brief Base building block for the implementation of an equation
 */

#ifndef QUICC_EQUATIONS_IEQUATION_HPP
#define QUICC_EQUATIONS_IEQUATION_HPP

// First include
//

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/MatrixOperationsInternal.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Equations/EquationParameters.hpp"
#include "QuICC/Equations/CouplingFeature.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Equations/EquationData.hpp"
#include "QuICC/Variables/VariableRequirement.hpp"
#include "QuICC/Simulation/SimulationBoundary.hpp"
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"
#include "QuICC/SpectralKernels/ISpectralKernel.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TransformConfigurators/TransformPath.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base building block for the implementation of an equation
    */
   class IEquation : public EquationData
   {
      public:
         /**
          * @brief Simple constructor
          *
          * \param spEqParams Shared equation parameters
          */
         explicit IEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IEquation();

         /**
          * @brief Access the shared resolution
          */
         virtual SharedResolution spRes() const = 0;

         /**
          * @brief Access the shared resolution
          */
         virtual const Resolution& res() const = 0;

         /**
          * @brief Initialise the equation
          */
         virtual void init(const SharedSimulationBoundary spBcIds);

         /**
          * @brief Initialise the nonlinear kernel
          */
         virtual void initNLKernel(const bool force = false);

         /**
          * @brief Get forward transform paths
          */
         virtual std::vector<Transform::TransformPath> forwardPaths();

         /**
          * @brief Get backward transform paths
          */
         virtual std::vector<Transform::TransformPath> backwardPaths() = 0;

         /**
          * @brief Generic model operator dispatcher to python scripts
          */
         virtual void buildModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id comp, const int matIdx, const std::size_t bcType) const; // = 0;

         /**
          * @brief Get nonlinear kernel
          */
         Physical::Kernel::SharedIPhysicalKernel spNLKernel() const;

         /**
          * @brief Initialise the spectral equation matrices
          */
         virtual void initSpectralMatrices() = 0;

         /**
          * @brief Implementation of the galerkin stencil dispatch to python scripts
          */
         void dispatchGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare = false) const;

         /**
          * @brief Set spectral constraint kernel
          */
         virtual void setConstraintKernel(FieldComponents::Spectral::Id compId, Spectral::Kernel::SharedISpectralKernel spKernel);

         /**
          * @brief Get spectral constraint kernel
          */
         Spectral::Kernel::SharedISpectralKernel spConstraintKernel(FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Initialize constraint kernels
          */
         virtual void initConstraintKernel();

         /**
          * @brief Set spectral source kernel
          */
         virtual void setSrcKernel(FieldComponents::Spectral::Id compId, Spectral::Kernel::SharedISpectralKernel spKernel);

         /**
          * @brief Get spectral source kernel
          */
         Spectral::Kernel::SharedISpectralKernel spSrcKernel(FieldComponents::Spectral::Id compId) const;

         /**
          * @brief Initialize source kernels
          */
         virtual void initSrcKernel();

      protected:
         /**
          * @brief Set the equation variable requirements
          */
         virtual void setRequirements() = 0;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() = 0;

         /**
          * @brief Set the default nonlinear components
          */
         virtual void setNLComponents() = 0;

         /**
          * @brief Initialise the spectral equation matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initSpectralMatricesComponent(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

         /**
          * @brief Implementation of the coupling definition to python scripts
          */
         void dispatchCoupling(FieldComponents::Spectral::Id comp, CouplingInformation::EquationTypeId eqType, const int iZero, const std::map<CouplingFeature,bool>& features, const Resolution& res);

         /**
          * @brief Implementation of model operator dispatcher to python scripts
          */
         void dispatchModelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, FieldComponents::Spectral::Id comp, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs) const;

         /**
          * @brief Implementation of the explicit matrix operator dispatch to python scripts
          */
         void dispatchExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId, const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs) const;

         /**
          * @brief Shared physical interaction kernel
          */
         Physical::Kernel::SharedIPhysicalKernel mspNLKernel;

         /**
          * @brief Shared spectral source kernel for each component
          */
          std::map<FieldComponents::Spectral::Id,Spectral::Kernel::SharedISpectralKernel> mSrcKernel;

         /**
          * @brief Shared spectral constraint kernel for each component
          */
          std::map<FieldComponents::Spectral::Id,Spectral::Kernel::SharedISpectralKernel> mConstraintKernel;

      private:

         /**
          * @brief Initialise the Galerkin stencisl
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initGalerkinStencils(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

         /**
          * @brief Initialise the quasi-inverse spectral  matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          */
         void initQIMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId);

         /**
          * @brief Initialise the explicit spectral  matrices for given component
          *
          * @param spBcIds List of boundary condition IDs
          * @param compId  Spectral component
          * @param opId    Type of explicit operator
          */
         void initExplicitMatrices(const SharedSimulationBoundary spBcIds, FieldComponents::Spectral::Id compId, const std::size_t opId);

         /**
          * @brief Set the galerkin stencil
          */
         virtual void setGalerkinStencil(FieldComponents::Spectral::Id compId, SparseMatrix &mat, const int matIdx) const; // = 0;

         /**
          * @brief Set the explicit matrix operator
          */
         virtual void setExplicitBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const std::size_t opId, const SpectralFieldId fieldId, const int matIdx) const; // = 0;

   };

   /// Typedef for a smart IEquation
   typedef std::shared_ptr<IEquation> SharedIEquation;

   /**
    * @brief Apply the quasi-inverse operator
    *
    * @param eq         Equation
    * @param compId     Equation field component ID
    * @param rField     Output field
    * @param start      Start index in linear storage
    * @param matIdx     System index
    * @param rhsStart   Start index in RHS data
    * @param rhs        RHS field data
    */
   template <typename TData> void applyQuasiInverse(const IEquation& eq, TData& rField, const int start, const int matIdx, const int rhsStart, const TData& rhs);
   template <> void applyQuasiInverse<DecoupledZMatrix>(const IEquation& eq, DecoupledZMatrix& rField, const int start, const int matIdx, const int rhsStart, const DecoupledZMatrix& rhs);

   /**
    * @brief Apply the galerkin stencil operator
    *
    * @param eq         Equation
    * @param compId     Equation field component ID
    * @param rField     Output field
    * @param start      Start index in linear storage
    * @param matIdx     System index
    * @param rhs        RHS field data
    */
   template <typename TData> void applyGalerkinStencil(const IEquation& eq, TData& rField, const int start, const int matIdx, const TData& rhs);
   template <> void applyGalerkinStencil<DecoupledZMatrix>(const IEquation& eq, DecoupledZMatrix& rField, const int start, const int matIdx, const DecoupledZMatrix& rhs);

   template <typename TData> inline void applyQuasiInverse(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& rField, const int start, const int matIdx, const int rhsStart, const TData& rhs)
   {
      if(eq.hasQID(compId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.quasiInverse<SparseMatrix>(compId, matIdx);

         // Get number of rows and cols
         int cols = rField.cols();
         int rhsRows = op->cols();

         Datatypes::internal::addMatrixProduct(rField, start, *op, rhs.block(rhsStart, 0, rhsRows, cols));

      } else if(eq.hasQIZ(compId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.quasiInverse<SparseMatrixZ>(compId, matIdx);

         // Get number of rows and cols
         int cols = rField.cols();
         int rhsRows = op->cols();

         Datatypes::internal::addMatrixProduct(rField, start, *op, rhs.block(rhsStart, 0, rhsRows, cols));
      }
   }

   template <> inline void applyQuasiInverse<DecoupledZMatrix>(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& rField, const int start, const int matIdx, const int rhsStart, const DecoupledZMatrix& rhs)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      if(eq.hasQID(compId))
      {
         // Create pointer to sparse operator
         const SparseMatrix * op = &eq.quasiInverse<SparseMatrix>(compId, matIdx);

         // Get number of rows and cols
         int cols = rField.real().cols();
         int rhsRows = op->cols();

         Datatypes::internal::addMatrixProduct(rField, start, *op, rhs.real().block(rhsStart, 0, rhsRows, cols), rhs.imag().block(rhsStart, 0, rhsRows, cols));

      } else if(eq.hasQIZ(compId))
      {
         // Create pointer to sparse operator
         const SparseMatrixZ * op = &eq.quasiInverse<SparseMatrixZ>(compId, matIdx);

         // Get number of rows and cols
         int cols = rField.real().cols();
         int rhsRows = op->cols();

         Datatypes::internal::addMatrixProduct(rField, start, *op, rhs.real().block(rhsStart, 0, rhsRows, cols), rhs.imag().block(rhsStart, 0, rhsRows, cols));
      }
   }

   template <typename TData> inline void applyGalerkinStencil(const IEquation& eq, FieldComponents::Spectral::Id compId, TData& rField, const int start, const int matIdx, const TData& rhs)
   {
      // Create pointer to sparse operator
      const SparseMatrix * op = &eq.galerkinStencil(compId, matIdx);

      Datatypes::internal::setMatrixProduct(rField, 0, *op, rhs.block(start, 0, op->cols(), rhs.cols()));
   }

   template <> inline void applyGalerkinStencil<DecoupledZMatrix>(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZMatrix& rField, const int start, const int matIdx, const DecoupledZMatrix& rhs)
   {
      assert(rField.real().rows() == rField.imag().rows());
      assert(rField.real().cols() == rField.imag().cols());

      // Create pointer to sparse operator
      const SparseMatrix * op = &eq.galerkinStencil(compId, matIdx);

      Datatypes::internal::setMatrixProduct(rField, 0, *op, rhs.real().block(start, 0, op->cols(), rhs.real().cols()), rhs.imag().block(start, 0, op->cols(), rhs.imag().cols()));
   }

   template <typename T, typename TOperator,typename TData> void computeExplicitTerm(const IEquation& eq, const std::size_t opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField, const int matIdx)
   {
      // Create pointer to sparse operator
      const TOperator * op = &eq.explicitOperator<TOperator>(opId, compId, fieldId, matIdx);

      if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         typename Eigen::Matrix<T,Eigen::Dynamic,1>  tmp(op->cols());
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            // Initialise storage to zero
            tmp.setZero();
            int l;
            int j_;
            int dimI = eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            int corrDim;
            if((eq.res().sim().ss().has(SpatialScheme::Feature::ShellGeometry) || eq.res().sim().ss().has(SpatialScheme::Feature::SphereGeometry)) &&
                  eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123) &&
                  eq.res().sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
            {
               corrDim = eq.res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
            }
            for(int j = 0; j < explicitField.slice(matIdx).cols(); j++)
            {
               j_ = eq.res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
               if(corrDim > 0)
               {
                  j_ -= corrDim;
               }
               for(int i = 0; i < explicitField.slice(matIdx).rows(); i++)
               {
                  // Compute correct position
                  l = j_ + i;

                  // Copy field value into storage
                  tmp(l) = explicitField.point(i,j,matIdx);
               }
            }
         #else
            int k = 0;
            for(int j = 0; j < explicitField.slice(matIdx).cols(); j++)
            {
               for(int i = 0; i < explicitField.slice(matIdx).rows(); i++)
               {
                  // Copy slice into flat array
                  tmp(k) = explicitField.point(i,j,matIdx);

                  // increase storage counter
                  k++;
               }
            }
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

         // Apply operator to field
         Datatypes::internal::addMatrixProduct(rSolverField, eqStart, *op, tmp);

      } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         // Apply operator to field
         Datatypes::internal::addMatrixProduct(rSolverField, eqStart, *op, explicitField.slice(matIdx));

      } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         // Get mode indexes
         ArrayI mode = eq.res().cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

         // Assert correct sizes
         assert(op->cols() == explicitField.slice(mode(0)).rows());

         // Apply operator to field
         Datatypes::internal::addMatrixProduct(rSolverField, eqStart, *op, explicitField.slice(mode(0)).col(mode(1)));

      } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         assert(matIdx == 0);

         /// \mhdBug very bad and slow implementation!
         typename Eigen::Matrix<T,Eigen::Dynamic,1>  tmp(op->cols());
         int l = 0, k_, j_, dimK, dimJ;

         switch(eq.res().sim().ss().dimension())
         {
            case 3:
               dimK = eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*eq.res().sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
               dimJ = eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               break;
            case 2:
               dimK = 1;
               dimJ = eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               break;
            case 1:
               dimK = 1;
               dimJ = 1;
               break;
         }

         for(int k = 0; k < eq.res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); k++)
         {
            k_ = eq.res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)*dimK;
            for(int j = 0; j < explicitField.slice(k).cols(); j++)
            {
               j_ = eq.res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
               for(int i = 0; i < explicitField.slice(k).rows(); i++)
               {
                  // Compute correct position
                  l = k_ + j_ + i;

                  // Copy slice into flat array
                  tmp(l) = explicitField.point(i,j,k);
               }
            }
         }

         // Apply operator to field
         Datatypes::internal::addMatrixProduct(rSolverField, eqStart, *op, tmp);
      }
   }

}
}

#endif // QUICC_EQUATIONS_IEQUATION_HPP
