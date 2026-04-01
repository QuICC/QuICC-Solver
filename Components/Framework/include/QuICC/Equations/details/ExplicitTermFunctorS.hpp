/**
 * @file ExplicitTermFunctorS.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORS_HPP
#define QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORS_HPP

// System includes
//

// Project includes
//
#include "Arithmetics/LinearAlgebra.hpp"
#include "Arithmetics/Utility.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/details/ExplicitTermFunctor.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <>
template <typename T, typename TOperator, typename TData>
void ExplicitTermFunctor<CouplingIndexType::SINGLE>::apply(
   TData& rSolverField, const TOperator& op, const int eqStart,
   const typename Framework::Selector::ScalarField<T>& explicitField)
{
   if constexpr ((std::is_same<T, MHDFloat>::value ||
                    std::is_same<T, MHDComplex>::value) &&
                 (std::is_same<TOperator, SparseMatrixZ>::value &&
                    std::is_same<TData, Matrix>::value))
   {}
   else
   {
      const auto& tRes = *res.cpu()->dim(Dimensions::Transform::SPECTRAL);
      const auto& sRes = res.sim();
      assert(matIdx == 0);

      /// \mhdBug very bad and slow implementation!
      typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp(op.cols(),
         1);
      int l = 0, k_, j_, dimK, dimJ;

      switch (sRes.ss().dimension())
      {
      case 3:
         dimK = sRes.dim(Dimensions::Simulation::SIM1D,
                   Dimensions::Space::SPECTRAL) *
                sRes.dim(Dimensions::Simulation::SIM3D,
                   Dimensions::Space::SPECTRAL);
         dimJ = sRes.dim(Dimensions::Simulation::SIM1D,
            Dimensions::Space::SPECTRAL);
         break;
      case 2:
         dimK = 1;
         dimJ = sRes.dim(Dimensions::Simulation::SIM1D,
            Dimensions::Space::SPECTRAL);
         break;
      case 1:
         dimK = 1;
         dimJ = 1;
         break;
      default:
         dimK = -1;
         dimJ = -1;
         throw std::logic_error("Spatial scheme has unknown dimension!");
      }

      for (int k = 0; k < tRes.template dim<Dimensions::Data::DAT3D>(); k++)
      {
         k_ = tRes.template idx<Dimensions::Data::DAT3D>(k) * dimK;
         const int cols = tRes.dim<Dimensions::Data::DAT2D>(matIdx);
         for (int j = 0; j < cols; j++)
         {
            j_ = tRes.template idx<Dimensions::Data::DAT2D>(j, k) * dimJ;
            const int usedRows =
               tRes.template dim<Dimensions::Data::DATB1D>(j, k);
            for (int i = 0; i < usedRows; i++)
            {
               // Compute correct position
               l = k_ + j_ + i;

               // Copy slice into flat array
               tmp(l, 0) = explicitField.point(i, j, k);
            }
         }
      }

      // Apply operator to field
      std::tuple<int, int, int, int> outBlk =
         std::make_tuple(eqStart, 0, op.rows(), Arithmetics::getCols(tmp));
      Arithmetics::computeAx<Arithmetics::Operation::Plus>(rSolverField, outBlk,
         op, tmp);
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTORS_HPP
