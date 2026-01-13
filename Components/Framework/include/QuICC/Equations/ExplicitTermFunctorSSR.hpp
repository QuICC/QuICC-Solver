/**
 * @file ExplicitTermFunctorSSR.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_EXPLICITTERMFUNCTORSSR_HPP
#define QUICC_EQUATIONS_EXPLICITTERMFUNCTORSSR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"
#include "QuICC/Equations/ExplicitTermFunctor.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "Arithmetics/Basic.hpp"

namespace QuICC {

namespace Equations {

template <>
   template <typename T, typename TOperator, typename TData>
      void ExplicitTermFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS>::compute(const std::size_t opId, TData& rSolverField, const int eqStart, SpectralFieldId fieldId, const typename Framework::Selector::ScalarField<T>& explicitField)
      {
         if constexpr((std::is_same<T,MHDFloat>::value || std::is_same<T, MHDComplex>::value ) && (std::is_same<TOperator, SparseMatrixZ>::value && std::is_same<TData, Matrix>::value))
         {
         } else
         {
            // Create pointer to sparse operator
            const TOperator * op = &eq->template explicitOperator<TOperator>(opId, compId, fieldId, matIdx);

            const auto& tRes = *eq->res().cpu()->dim(Dimensions::Transform::SPECTRAL);
            typename Eigen::Matrix<T,Eigen::Dynamic,1>  tmp(op->cols());
#if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            const auto& sRes = eq->res().sim();
            // Initialise storage to zero
            tmp.setZero();
            int l;
            int j_;
            int dimI = sRes.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
            int corrDim;
            if((sRes.ss().has(SpatialScheme::Feature::ShellGeometry) || sRes.ss().has(SpatialScheme::Feature::SphereGeometry)) &&
                  sRes.ss().has(SpatialScheme::Feature::SpectralOrdering123) &&
                  sRes.ss().has(SpatialScheme::Feature::SpectralMatrix2D))
            {
               corrDim = tRes.template idx<Dimensions::Data::DAT3D>(matIdx)*dimI;
            }
            const int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
            for(int j = 0; j < cols; j++)
            {
               j_ = tRes.template idx<Dimensions::Data::DAT2D>(j,matIdx)*dimI;
               if(corrDim > 0)
               {
                  j_ -= corrDim;
               }
               int rows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);
               for(int i = 0; i < rows; i++)
               {
                  // Compute correct position
                  l = j_ + i;

                  // Copy field value into storage
                  tmp(l) = explicitField.point(i,j,matIdx);
               }
            }
#else
            int k = 0;
            const int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
            for(int j = 0; j < cols; j++)
            {
               // Effective rows in case of non-uniform truncation
               int usedRows = tRes.dim<Dimensions::Data::DATB1D>(j, matIdx);

               for(int i = 0; i < usedRows; i++)
               {
                  // Copy slice into flat array
                  tmp(k) = explicitField.point(i,j,matIdx);

                  // increase storage counter
                  k++;
               }
            }
#endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

            // Apply operator to field
            Arithmetics::addMatrixProduct(rSolverField, eqStart, *op, tmp);
         }
      }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_EXPLICITTERMFUNCTORSSR_HPP
