/**
 * @file ExplicitTermFunctor.hpp
 * @brief ExplicitTerm calculation implementation
 */

#ifndef QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTOR_HPP
#define QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class ExplicitTermFunctor
{
public:
   /**
    * @brief ctor
    */
   ExplicitTermFunctor(const Resolution& res, const CouplingInformation& cinfo,
      const int matIdx, const bool useShift, const bool shiftTop);

   /**
    * @brief deleted default ctor
    */
   ExplicitTermFunctor() = delete;

   /**
    * @brief dtor
    */
   ~ExplicitTermFunctor() = default;

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param rSolverField  Solver field values
    * @param eqStart       Start index for the equation field
    * @param explicitField Explicit linear field values
    */
   template <typename T, typename TOperator, typename TData>
   void apply(TData& rSolverField, const TOperator& op, const int eqStart,
      const Framework::Selector::ScalarField<T>& explicitField);

private:
   /**
    * @brief Initialization
    */
   void init(const bool shiftTop);

   /**
    * @brief Resolution
    */
   const Resolution& res;

   /**
    * @brief Coupling information
    */
   const CouplingInformation& cinfo;

   /**
    * @brief Matrix index
    */
   const int matIdx;

   /**
    * @brief Starting row
    */
   int zeroRow;

   /**
    * @brief Starting column
    */
   int zeroCol;

   /**
    * @brief Shift max row
    */
   int shiftMaxRow;

   /**
    * @brief Shift max col
    */
   int shiftMaxCol;
};

template <CouplingIndexType IndexType>
void ExplicitTermFunctor<IndexType>::init(const bool shiftTop)
{
   zeroRow = cinfo.galerkinShift(matIdx, 0);
   if (res.sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
   {
      zeroCol = cinfo.galerkinShift(matIdx, 2);
   }
   else
   {
      zeroCol = cinfo.galerkinShift(matIdx, 1);
   }
   shiftMaxRow = 0;
   shiftMaxCol = 0;

   if(!shiftTop)
   {
      std::swap(shiftMaxRow, zeroRow);
      std::swap(shiftMaxCol, zeroCol);
   }
}

template <CouplingIndexType IndexType>
ExplicitTermFunctor<IndexType>::ExplicitTermFunctor(const Resolution& res, const CouplingInformation& cinfo,
   const int matIdx, const bool useShift, const bool shiftTop) :
    res(res), cinfo(cinfo), matIdx(matIdx), zeroRow(0), zeroCol(0), shiftMaxRow(0), shiftMaxCol(0)
{
   if (useShift)
   {
      this->init(shiftTop);
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_EXPLICITTERMFUNCTOR_HPP
