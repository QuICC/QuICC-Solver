/**
 * @file CopyUnknownFunctor.hpp
 * @brief Implementation of the CopyUnknown functors
 */

#ifndef QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTOR_HPP
#define QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/Equations/IFieldEquation.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class CopyUnknownFunctor
{
public:
   /**
    * @brief ctor
    */
   CopyUnknownFunctor(const IFieldEquation& eq,
      FieldComponents::Spectral::Id compId, const int matIdx,
      const bool useShift, const bool shiftTop);

   /**
    * @brief deleted default ctor
    */
   CopyUnknownFunctor() = delete;

   /**
    * @brief dtor
    */
   ~CopyUnknownFunctor() = default;

   /**
    * @brief Apply functor
    */
   template <bool IsSet, typename TData, typename TField>
   void apply(const TField& field, TData& storage, const int start);

private:
   /**
    * @brief Initialization
    */
   void init(const bool shiftTop);

   /**
    * @brief Reference to equation
    */
   const IFieldEquation* eq;

   /**
    * @brief Field component ID
    */
   FieldComponents::Spectral::Id compId;

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
void CopyUnknownFunctor<IndexType>::init(const bool shiftTop)
{
   const auto& info = eq->couplingInfo(compId);
   zeroRow = info.galerkinShift(matIdx, 0);
   if (eq->res().sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
   {
      zeroCol = info.galerkinShift(matIdx, 2);
   }
   else
   {
      zeroCol = info.galerkinShift(matIdx, 1);
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
CopyUnknownFunctor<IndexType>::CopyUnknownFunctor(const IFieldEquation& eq,
   FieldComponents::Spectral::Id compId, const int matIdx,
   const bool useShift, const bool shiftTop) :
    eq(&eq), compId(compId), matIdx(matIdx), zeroRow(0), zeroCol(0), shiftMaxRow(0), shiftMaxCol(0)
{
   if (useShift)
   {
      this->init(shiftTop);
   }
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTOR_HPP
