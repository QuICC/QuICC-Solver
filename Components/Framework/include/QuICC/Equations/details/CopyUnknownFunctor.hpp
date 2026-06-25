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
#include "QuICC/Equations/CouplingFeature.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class CopyUnknownFunctor
{
public:
   /**
    * @brief ctor
    */
   CopyUnknownFunctor(const Resolution& res, const Equations::CouplingInformation& cinfo,
      const int matIdx,
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
    * @brief Resolution
    */
   const Resolution& res;

   /**
    * @brief Coupling information
    */
   const Equations::CouplingInformation& cinfo;

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
CopyUnknownFunctor<IndexType>::CopyUnknownFunctor(const Resolution& res, const Equations::CouplingInformation& cinfo,
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

#endif // QUICC_EQUATIONS_DETAILS_COPYUNKNOWNFUNCTOR_HPP
