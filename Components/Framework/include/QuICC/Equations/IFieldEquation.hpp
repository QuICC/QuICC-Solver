/**
 * @file IFieldEquation.hpp
 * @brief Base building block for the implementation of an equation
 */

#ifndef QUICC_EQUATIONS_IFIELDEQUATION_HPP
#define QUICC_EQUATIONS_IFIELDEQUATION_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IEquation.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/DecoupledComplexInternal.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Base building block for the implementation of an equation
    */
   class IFieldEquation : public IEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          */
         explicit IFieldEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple constructor
          *
          * @param spEqParams Equation parameters
          * @param spScheme   Spatial scheme
          * @param spBackend  Model backend
          * @param spOptions  Additional options
          */
         explicit IFieldEquation(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<EquationOptions> spOptions);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IFieldEquation() = default;

         /**
          * @brief Apply generic constraint on spectral data
          *
          * @param compId  ID of the spectral component
          * @param timeId  Timing of the constraint
          * @return contraint changed solution?
          */
         virtual bool applyConstraint(FieldComponents::Spectral::Id compId, const std::size_t timeId);

         /**
          * @brief Generic source term implementation
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual MHDVariant sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

         /**
          * @brief Generic boundary value implementation
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual MHDVariant boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

      protected:
         /**
          * @brief Templated passthrough update the stored value with the solver solution
          */
         virtual MHDVariant updateStoredSolution(const MHDVariant newData, FieldComponents::Spectral::Id compId, const int i, const int j, const int k);

         /**
          * @brief Transfer solver solution to equation unknown
          *
          * @param field   Scalar or vector field
          * @param compId  Component ID
          * @param storage Solver solution
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         template <typename TData, typename TField> void storeSolutionImpl(TField& field, FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start);
      private:

   };

   /// Typedef for a smart IFieldEquation
   typedef std::shared_ptr<IFieldEquation> SharedIFieldEquation;

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param eq            Equation
    * @param compId        Equation field component ID
    * @param rSolverField  Solver field values
    * @param eqStart       Start index for the equation field
    * @param fieldId       Physical field ID
    * @param explicitField Explicit linear field values
    * @param matIdx        System index
    */
   template <typename T, typename TData>
      void addExplicitTerm(const IFieldEquation& eq, const std::size_t opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField, const int matIdx);
   template <typename T, typename TOperator,typename TData> void computeExplicitTerm(const IFieldEquation& eq, const std::size_t opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const Framework::Selector::ScalarField<T>& explicitField, const int matIdx);

   /**
    * @brief Copy unknown spectral values to solver
    *
    * @param eq         Equation to work on
    * @param compId     Component ID
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    * @param useShift   Use galerkin shifts
    * @param isSet      Arithmetic operation is set
    */
   template <typename TData, typename TField> void copyUnknown(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet);

   /**
    * @brief Add source term
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField> void addSource(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Set boundary value
    *
    * @param eq      Equation to work on
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField> void setBoundaryValue(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);

   /**
    * @brief Set nonlinear spectral values to zero
    *
    * @param compId  Component ID
    * @param storage Storage for the equation values
    * @param matIdx  Index of the given data
    * @param start   Start index for the storage
    */
   template <typename TData, typename TField> void setZeroNonlinear(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start);


   template <typename T, typename TData>
      void addExplicitTerm(const IFieldEquation& eq, const std::size_t opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField, const int matIdx)
   {
      // Compute with complex linear operator
      if(eq.hasExplicitZTerm(opId, compId, fieldId))
      {
         computeExplicitTerm<T,SparseMatrixZ>(eq, opId, compId, rSolverField,  eqStart, fieldId, explicitField, matIdx);
      }

      // Compute with real linear operator
      if(eq.hasExplicitDTerm(opId, compId, fieldId))
      {
         computeExplicitTerm<T,SparseMatrix>(eq, opId, compId, rSolverField,  eqStart, fieldId, explicitField, matIdx);
      }
   }

   template <typename T, typename TOperator, typename TData>
      void computeExplicitTerm(const IFieldEquation& eq, const std::size_t opId, FieldComponents::Spectral::Id compId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField, const int matIdx)
   {
      if constexpr((std::is_same<T,MHDFloat>::value || std::is_same<T, MHDComplex>::value ) && (std::is_same<TOperator, SparseMatrixZ>::value && std::is_same<TData, Matrix>::value))
      {
      } else
      {
         // Create pointer to sparse operator
         const TOperator * op = &eq.template explicitOperator<TOperator>(opId, compId, fieldId, matIdx);

         const auto& tRes = *eq.res().cpu()->dim(Dimensions::Transform::SPECTRAL);
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
                  corrDim = tRes.template idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
               }
               for(int j = 0; j < explicitField.slice(matIdx).cols(); j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
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
            ArrayI mode = tRes.mode(matIdx);

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
               default:
                  dimK = -1;
                  dimJ = -1;
                  throw  std::logic_error("Spatial scheme has unknown dimension!");
            }

            for(int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = tRes.template idx<Dimensions::Data::DAT3D>(k)*dimK;
               for(int j = 0; j < explicitField.slice(k).cols(); j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
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

   template <typename TData, typename TField> void copyUnknown(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet)
   {
      int zeroRow = 0;
      int zeroCol = 0;
      if(useShift)
      {
         zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
         if(eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
         {
            zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
         } else
         {
            zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
         }
      }

      const auto& tRes = *eq.res().cpu()->dim(Dimensions::Transform::SPECTRAL);
      // matIdx is the index of the slowest varying direction with a single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         int rows = field.comp(compId).slice(matIdx).rows();
         int cols = field.comp(compId).slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);

         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            // Add source data
            int l;
            int j_;
            int dimI = eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            int corrDim;
            if((eq.res().sim().ss().has(SpatialScheme::Feature::ShellGeometry) || eq.res().sim().ss().has(SpatialScheme::Feature::SphereGeometry)) &&
                  eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123) &&
                  eq.res().sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
            {
               corrDim = tRes.template idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
            }
            if(isSet)
            {
               ///\mhdBug This is overkill
               // Set storage to zero
               setZeroNonlinear(eq, field, compId, storage, matIdx, start);

               for(int j = zeroCol; j < cols; j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  if(corrDim > 0)
                  {
                     j_ -= corrDim;
                  }
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Copy field value into storage
                     Datatypes::internal::setScalar(storage, l, field.comp(compId).point(i,j,matIdx));
                  }
               }
            } else
            {
               for(int j = zeroCol; j < cols; j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  if(corrDim > 0)
                  {
                     j_ -= corrDim;
                  }
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Copy field value into storage
                     Datatypes::internal::addScalar(storage, l, field.comp(compId).point(i,j,matIdx));
                  }
               }
            }
         #else
            // Copy data
            int k = start;
            if(isSet)
            {
               for(int j = zeroCol; j < cols; j++)
               {
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Copy field value into storage
                     Datatypes::internal::setScalar(storage, k, field.comp(compId).point(i,j,matIdx));

                     // increase storage counter
                     k++;
                  }
               }
            } else
            {
               for(int j = zeroCol; j < cols; j++)
               {
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Copy field value into storage
                     Datatypes::internal::addScalar(storage, k, field.comp(compId).point(i,j,matIdx));

                     // increase storage counter
                     k++;
                  }
               }
            }
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

      // matIdx is the index of the slowest varying direction with multiple RHS
      } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         int rows = field.comp(compId).slice(matIdx).rows();
         int cols = field.comp(compId).slice(matIdx).cols();

         //Safety assertion
         assert(start >= 0);

         // Copy data
         if(isSet)
         {
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Copy field value into storage
                  Datatypes::internal::setScalar(storage, i - zeroRow + start, j - zeroCol, field.comp(compId).point(i,j,matIdx));
               }
            }
         } else
         {
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Copy field value into storage
                  Datatypes::internal::addScalar(storage, i - zeroRow + start, j - zeroCol, field.comp(compId).point(i,j,matIdx));
               }
            }
         }

      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         //Safety assertion
         assert(start >= 0);

         // Get mode indexes
         ArrayI mode = tRes.mode(matIdx);
         int rows = field.comp(compId).slice(mode(0)).rows();

         bool isUsed;
         if(eq.res().sim().ss().has(SpatialScheme::Feature::FourierIndex23))
         {
            // Filter out complex conjugate modes to be safe
            isUsed = !(mode(3) == 0 && mode(2) > eq.res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2);
         } else
         {
            isUsed = true;
         }

         if(isUsed)
         {
            // Copy data
            int k = start;
            if(isSet)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Copy field value into storage
                  Datatypes::internal::setScalar(storage, k, field.comp(compId).point(i,mode(1),mode(0)));

                  // increase storage counter
                  k++;
               }
            } else
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Copy field value into storage
                  Datatypes::internal::addScalar(storage, k, field.comp(compId).point(i,mode(1),mode(0)));

                  // increase storage counter
                  k++;
               }
            }

         } else
         {
            setZeroNonlinear(eq, field, compId, storage, matIdx, start);
         }

      // There is a single matrix
      } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         assert(matIdx == 0);

         //Safety assertion
         assert(start >= 0);

         // Copy data
         int l, k_, j_, dimK, dimJ;

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
            default:
               dimK = -1;
               dimJ = -1;
               throw  std::logic_error("Spatial scheme has unknown dimension!");
         }

         if(isSet)
         {
            for(int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = tRes.template idx<Dimensions::Data::DAT3D>(k)*dimK;
               for(int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
                  for(int i = 0; i < eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
                  {
                     // Compute correct position
                     l = start + k_ + j_ + i;

                     // Copy field value into storage
                     Datatypes::internal::setScalar(storage, l, field.comp(compId).point(i,j,k));
                  }
               }
            }
         } else
         {
            for(int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = tRes.template idx<Dimensions::Data::DAT3D>(k)*dimK;
               for(int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
                  for(int i = 0; i < eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
                  {
                     // Compute correct position
                     l = start + k_ + j_ + i;

                     // Copy field value into storage
                     Datatypes::internal::addScalar(storage, l, field.comp(compId).point(i,j,k));
                  }
               }
            }
         }
      }
   }

   template <typename TData, typename TField> void addSource(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Add source term if required
      if(eq.couplingInfo(compId).hasSource())
      {
         const auto& tRes = *eq.res().cpu()->dim(Dimensions::Transform::SPECTRAL); 
         // matIdx is the index of the slowest varying direction with a single RHS
         if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
         {
            int rows = field.comp(compId).slice(matIdx).rows();
            int cols = field.comp(compId).slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            int zeroCol;
            if(eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
            } else
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            }

            //Safety assertion
            assert(start >= 0);

            #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
               // Add source data
               int l;
               int j_;
               int dimI = eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               int corrDim;
               if((eq.res().sim().ss().has(SpatialScheme::Feature::ShellGeometry) || eq.res().sim().ss().has(SpatialScheme::Feature::SphereGeometry)) &&
                     eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123) &&
                     eq.res().sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
               {
                  corrDim = tRes.template idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
               }
               for(int j = zeroCol; j < cols; j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  if(corrDim > 0)
                  {
                     j_ -= corrDim;
                  }
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Add source term
                     Datatypes::internal::addScalar(storage, l, eq.sourceTerm(compId, i, j, matIdx));
                  }
               }
            #else
               // Add source term
               int k = start;
               for(int j = zeroCol; j < cols; j++)
               {
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Add source term
                     Datatypes::internal::addScalar(storage, k, eq.sourceTerm(compId, i, j, matIdx));

                     // increase storage counter
                     k++;
                  }
               }
            #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

         // matIdx is the index of the slowest varying direction with multiple RHS
         } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
         {
            int rows = field.comp(compId).slice(matIdx).rows();
            int cols = field.comp(compId).slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            int zeroCol;
            if(eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
            } else
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            }

            //Safety assertion
            assert(start >= 0);

            // Copy data
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Add source term
                  Datatypes::internal::addScalar(storage, i - zeroRow + start, j - zeroCol, eq.sourceTerm(compId, i, j, matIdx));
               }
            }

         // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
         {
            //Safety assertion
            assert(start >= 0);

            // Get mode indexes
            ArrayI mode = tRes.mode(matIdx);
            int rows = field.comp(compId).slice(mode(0)).rows();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);

            // Copy data
            int k = start;
            for(int i = zeroRow; i < rows; i++)
            {
               // Add source term
               Datatypes::internal::addScalar(storage, k, eq.sourceTerm(compId, i, mode(1), mode(0)));

               // increase storage counter
               k++;
            }

         // There is a single matrix
         } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
         {
            assert(matIdx == 0);

            //int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            //int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            //int zeroBlock = eq.couplingInfo(compId).galerkinShift(matIdx,2);

            //Safety assertion
            assert(start >= 0);

            // Add source term
            int l, k_, j_, dimK, dimJ;

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
               default:
                  dimK = -1;
                  dimJ = -1;
                  throw  std::logic_error("Spatial scheme has unknown dimension!");
            }

            for(int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = tRes.template idx<Dimensions::Data::DAT3D>(k)*dimK;
               for(int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
                  for(int i = 0; i < eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
                  {
                     // Compute correct position
                     l = start + k_ + j_ + i;

                     // Add source term
                     Datatypes::internal::addScalar(storage, l, eq.sourceTerm(compId, i, j, k));
                  }
               }
            }
         }
      }
   }

   template <typename TData, typename TField> void setBoundaryValue(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      // Set boundary value if required
      if(eq.couplingInfo(compId).hasBoundaryValue())
      {
         if(eq.couplingInfo(compId).isGalerkin())
         {
            throw std::logic_error("Galerkin expansion cannot have a nonzero boundary value!");
         }

         const auto& tRes = *eq.res().cpu()->dim(Dimensions::Transform::SPECTRAL);
         // matIdx is the index of the slowest varying direction with a single RHS
         if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
         {
            int rows = field.comp(compId).slice(matIdx).rows();
            int cols = field.comp(compId).slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            int zeroCol;
            if(eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
            } else
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            }

            //Safety assertion
            assert(start >= 0);

            #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
               // Set boundary value
               int l;
               int j_;
               int dimI = eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               int corrDim;
               if((eq.res().sim().ss().has(SpatialScheme::Feature::ShellGeometry) || eq.res().sim().ss().has(SpatialScheme::Feature::SphereGeometry)) &&
                     eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123) &&
                     eq.res().sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
               {
                  corrDim = tRes.template idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
               }
               for(int j = zeroCol; j < cols; j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
                  if(corrDim > 0)
                  {
                     j_ -= corrDim;
                  }
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Compute correct position
                     l = start + j_ + i;

                     // Add source term
                     Datatypes::internal::setScalar(storage, l, eq.boundaryValue(compId, i, j, matIdx));
                  }
               }
            #else
               // Set boundary value
               int k = start;
               for(int j = zeroCol; j < cols; j++)
               {
                  for(int i = zeroRow; i < rows; i++)
                  {
                     // Add source term
                     Datatypes::internal::setScalar(storage, k, eq.boundaryValue(compId, i, j, matIdx));

                     // increase storage counter
                     k++;
                  }
               }
            #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

         // matIdx is the index of the slowest varying direction with multiple RHS
         } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
         {
            int rows = field.comp(compId).slice(matIdx).rows();
            int cols = field.comp(compId).slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            int zeroCol;
            if(eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
            } else
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            }

            //Safety assertion
            assert(start >= 0);

            // Copy data
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Add source term
                  Datatypes::internal::setScalar(storage, i - zeroRow + start, j - zeroCol, eq.boundaryValue(compId, i, j, matIdx));
               }
            }

         // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
         } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
         {
            //Safety assertion
            assert(start >= 0);

            // Get mode indexes
            ArrayI mode = tRes.mode(matIdx);
            int rows = field.comp(compId).slice(mode(0)).rows();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);

            // Copy data
            int k = start;
            for(int i = zeroRow; i < rows; i++)
            {
               // Add source term
               Datatypes::internal::setScalar(storage, k, eq.boundaryValue(compId, i, mode(1), mode(0)));

               // increase storage counter
               k++;
            }

         // There is a single matrix
         } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
         {
            assert(matIdx == 0);

            //int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            //int zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            //int zeroBlock = eq.couplingInfo(compId).galerkinShift(matIdx,2);

            //Safety assertion
            assert(start >= 0);

            // Add source term
            int l, k_, j_, dimK, dimJ;

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
               default:
                  dimK = -1;
                  dimJ = -1;
                  throw  std::logic_error("Spatial scheme has unknown dimension!");
            }

            for(int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = tRes.template idx<Dimensions::Data::DAT3D>(k)*dimK;
               for(int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
                  for(int i = 0; i < eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
                  {
                     // Compute correct position
                     l = start + k_ + j_ + i;

                     // Add source term
                     Datatypes::internal::setScalar(storage, l, eq.boundaryValue(compId, i, j, k));
                  }
               }
            }
         }
      }
   }

   template <typename TData, typename TField> void setZeroNonlinear(const IFieldEquation& eq, const TField& field, FieldComponents::Spectral::Id compId, TData& storage, const int matIdx, const int start)
   {
      const auto& tRes = *eq.res().cpu()->dim(Dimensions::Transform::SPECTRAL);
      // matIdx is the index of the slowest varying direction with a single RHS
      if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         //Safety assertion
         assert(start >= 0);

         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            for(int k = 0; k < eq.couplingInfo(compId).galerkinN(matIdx); ++k)
            {
               // Set field to zero
               Datatypes::internal::setScalar(storage, k + start, typename TData::Scalar(0.0));
            }
         #else
            int rows = field.comp(compId).slice(matIdx).rows();
            int cols = field.comp(compId).slice(matIdx).cols();
            int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
            int zeroCol;
            if(eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
            } else
            {
               zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
            }

            // Set data to zero
            int k = start;
            for(int j = zeroCol; j < cols; j++)
            {
               for(int i = zeroRow; i < rows; i++)
               {
                  // Set field to zero
                  Datatypes::internal::setScalar(storage, k, typename TData::Scalar(0.0));

                  // increase storage counter
                  k++;
               }
            }
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

      // matIdx is the index of the slowest varying direction with multiple RHS
      } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         int rows = field.comp(compId).slice(matIdx).rows();
         int cols = field.comp(compId).slice(matIdx).cols();
         int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);
         int zeroCol;
         if(eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
         {
            zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,2);
         } else
         {
            zeroCol = eq.couplingInfo(compId).galerkinShift(matIdx,1);
         }

         //Safety assertion
         assert(start >= 0);

         // Set data to zero
         for(int j = zeroCol; j < cols; j++)
         {
            for(int i = zeroRow; i < rows; i++)
            {
               // Set field to zero
               Datatypes::internal::setScalar(storage, i - zeroRow + start, j - zeroCol, typename TData::Scalar(0.0));
            }
         }

      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         //Safety assertion
         assert(start >= 0);

         // Get mode indexes
         ArrayI mode = tRes.mode(matIdx);
         int rows = field.comp(compId).slice(mode(0)).rows();
         int zeroRow = eq.couplingInfo(compId).galerkinShift(matIdx,0);

         // Set data to zero
         int k = start;
         for(int i = zeroRow; i < rows; i++)
         {
            // Set field to zero
            Datatypes::internal::setScalar(storage, k, typename TData::Scalar(0.0));

            // increase storage counter
            k++;
         }

      } else if(eq.couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         //Safety assertion
         assert(matIdx == 0);
         assert(start >= 0);

         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            for(int k = 0; k < eq.couplingInfo(compId).galerkinN(matIdx); ++k)
            {
               // Set field to zero
               Datatypes::internal::setScalar(storage, k + start, typename TData::Scalar(0.0));
            }
         #else
            // Set data to zero
            int l, k_, j_, dimK, dimJ;

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
               default:
                  dimK = -1;
                  dimJ = -1;
                  throw  std::logic_error("Spatial scheme has unknown dimension!");
            }

            for(int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = tRes.template idx<Dimensions::Data::DAT3D>(k)*dimK;
               for(int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
                  for(int i = 0; i < eq.res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
                  {
                     // Compute correct position
                     l = start + k_ + j_ + i;

                     // Set field to zero
                     Datatypes::internal::setScalar(storage, l, typename TData::Scalar(0.0));
                  }
               }
            }
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE
      }
   }

   inline MHDVariant IFieldEquation::updateStoredSolution(const MHDVariant newData, FieldComponents::Spectral::Id, const int, const int, const int)
   {
      return newData;
   }

   template<typename TData, typename TField>
      void IFieldEquation::storeSolutionImpl(TField& field, FieldComponents::Spectral::Id compId, const TData& storage, const int matIdx, const int start)
   {
      const TData * solution;
      TData tmp;
      int solStart;
      if(this->couplingInfo(compId).isGalerkin())
      {
         // Temporary storage is required
         tmp = TData(this->couplingInfo(compId).tauN(matIdx), this->couplingInfo(compId).rhsCols(matIdx));

         // Apply Galerkin stencil
         applyGalerkinStencil(*this, compId, tmp, start, matIdx, storage);

         solStart = 0;
         solution = &tmp;

      } else
      {
         solStart = start;
         solution = &storage;
      }

      const auto& tRes = *this->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
      if(this->couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         int rows = field.comp(compId).slice(matIdx).rows();
         int cols = field.comp(compId).slice(matIdx).cols();

         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            // Add source data
            int l;
            int j_;
            int dimI = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            int corrDim;
            if((eq.res().sim().ss().has(SpatialScheme::Feature::ShellGeometry) || eq.res().sim().ss().has(SpatialScheme::Feature::SphereGeometry)) &&
                  eq.res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering123) &&
                  eq.res().sim().ss().has(SpatialScheme::Feature::SpectralMatrix2D))
            {
               corrDim = tRes.template idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
            }
            for(int j = 0; j < cols; j++)
            {
               j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
               if(corrDim > 0)
               {
                  j_ -= corrDim;
               }
               for(int i = 0; i < rows; i++)
               {
                  // Compute correct position
                  l = start + j_ + i;

                  // Copy timestep output into field
                  MHDVariant dataPoint = Datatypes::internal::getScalar(*solution, l);
                  dataPoint = this->updateStoredSolution(dataPoint, compId, i, j, matIdx);
                  field.rComp(compId).setPoint(dataPoint,i,j,matIdx);
               }
            }
         #else
            // Copy data
            int k = solStart;
            for(int j = 0; j < cols; j++)
            {
               for(int i = 0; i < rows; i++)
               {
                  // Copy timestep output into field
                  MHDVariant dataPoint = Datatypes::internal::getScalar(*solution, k);
                  dataPoint = this->updateStoredSolution(dataPoint, compId, i, j, matIdx);
                  field.rComp(compId).setPoint(dataPoint,i,j,matIdx);

                  // increase linear storage counter
                  k++;
               }
            }
         #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

      } else if(this->couplingInfo(compId).indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         int rows = field.comp(compId).slice(matIdx).rows();
         int cols = field.comp(compId).slice(matIdx).cols();

         // Copy data
         for(int j = 0; j < cols; j++)
         {
            for(int i = 0; i < rows; i++)
            {
               // Copy timestep output into field
               MHDVariant dataPoint = Datatypes::internal::getScalar(*solution, i + solStart,j);
               dataPoint = this->updateStoredSolution(dataPoint, compId, i, j, matIdx);
               field.rComp(compId).setPoint(dataPoint,i,j,matIdx);
            }
         }

      } else if(this->couplingInfo(compId).indexType() == CouplingIndexType::MODE)
      {
         // Get mode indexes
         ArrayI mode = tRes.mode(matIdx);
         int rows = field.comp(compId).slice(mode(0)).rows();

         // Copy data
         int k = solStart;
         for(int i = 0; i < rows; i++)
         {
            // Copy timestep output into field
            MHDVariant dataPoint = Datatypes::internal::getScalar(*solution, k);
            dataPoint = this->updateStoredSolution(dataPoint, compId, i, mode(1), mode(0));
            field.rComp(compId).setPoint(dataPoint,i,mode(1),mode(0));

            // increase linear storage counter
            k++;
         }

      } else if(this->couplingInfo(compId).indexType() == CouplingIndexType::SINGLE)
      {
         assert(matIdx == 0);

         // Copy data
         int l, k_, j_, dimK, dimJ;

         switch(this->res().sim().ss().dimension())
         {
            case 3:
               dimK = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*this->res().sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
               dimJ = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               break;
            case 2:
               dimK = 1;
               dimJ = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
               break;
            case 1:
               dimK = 1;
               dimJ = 1;
               break;
            default:
               dimK = -1;
               dimJ = -1;
               throw  std::logic_error("Spatial scheme has unknown dimension!");
         }

         for(int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
         {
            k_ = tRes.template idx<Dimensions::Data::DAT3D>(k)*dimK;
            for(int j = 0; j < tRes.template dim<Dimensions::Data::DAT2D>(k); j++)
            {
               j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
               for(int i = 0; i < this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
               {
                  // Compute correct position
                  l = solStart + k_ + j_ + i;

                  // Copy timestep output into field
                  MHDVariant dataPoint = Datatypes::internal::getScalar(*solution, l);
                  dataPoint = this->updateStoredSolution(dataPoint, compId, i, j, k);
                  field.rComp(compId).setPoint(dataPoint,i,j,k);
               }
            }
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_IFIELDEQUATION_HPP
