/**
 * @file CopyUnknown.hpp
 * @brief CopyUnknown function
 */

#ifndef QUICC_EQUATIONS_COPYUNKNOWN_HPP
#define QUICC_EQUATIONS_COPYUNKNOWN_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "Arithmetics/Basic.hpp"
#include "QuICC/Equations/details/CopyUnknownFunctor.hpp"
#include "QuICC/Equations/details/CopyUnknownFunctorSMR.hpp"
#include "QuICC/Equations/details/CopyUnknownFunctorSSR.hpp"
#include "QuICC/Equations/details/CopyUnknownFunctorS.hpp"
#include "QuICC/Equations/details/CopyUnknownFunctorM.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Copy unknown spectral values to solver
    *
    * @param storage    Storage for the equation values
    * @param matIdx     Index of the given data
    * @param start      Start index for the storage
    * @param useShift   Use galerkin shifts
    * @param isSet      Arithmetic operation is set
    * @param shiftTop   shift top?
    */
   template <typename TData, typename TField> void copyUnknown(const Resolution& res, const CouplingInformation& cinfo, const TField& field, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet, const bool shiftTop);

   template <typename TData, typename TField> void copyUnknown(const Resolution& res, const CouplingInformation& cinfo, const TField& field, TData& storage, const int matIdx, const int start, const bool useShift, const bool isSet, const bool shiftTop)
   {
      // matIdx is the index of the slowest varying direction with a single RHS
      if(cinfo.indexType() == CouplingIndexType::SLOWEST_SINGLE_RHS)
      {
         details::CopyUnknownFunctor<CouplingIndexType::SLOWEST_SINGLE_RHS> func(res, cinfo, matIdx, useShift, shiftTop);

         // Copy data
         if(isSet)
         {
            func.apply<true>(field, storage, start);
         }
         else
         {
            func.apply<false>(field, storage, start);
         }
      }
      // matIdx is the index of the slowest varying direction with multiple RHS
      else if(cinfo.indexType() == CouplingIndexType::SLOWEST_MULTI_RHS)
      {
         details::CopyUnknownFunctor<CouplingIndexType::SLOWEST_MULTI_RHS> func(res, cinfo, matIdx, useShift, shiftTop);

         // Copy data
         if(isSet)
         {
            func.apply<true>(field, storage, start);
         }
         else
         {
            func.apply<false>(field, storage, start);
         }
      }
      // matIdx is the index of a 2D mode, conversion to the two (k,m) mode indexes required
      else if(cinfo.indexType() == CouplingIndexType::MODE)
      {
         details::CopyUnknownFunctor<CouplingIndexType::MODE> func(res, cinfo, matIdx, useShift, shiftTop);

         // Copy data
         if(isSet)
         {
            func.apply<true>(field, storage, start);
         }
         else
         {
            func.apply<false>(field, storage, start);
         }
      }
      // There is a single matrix
      else if(cinfo.indexType() == CouplingIndexType::SINGLE)
      {
         details::CopyUnknownFunctor<CouplingIndexType::SINGLE> func(res, cinfo, matIdx, useShift, shiftTop);

         // Copy data
         if(isSet)
         {
            func.apply<true>(field, storage, start);
         }
         else
         {
            func.apply<false>(field, storage, start);
         }
      }
   }

} // Equations
} // QuICC

#endif // QUICC_EQUATIONS_COPYUNKNOWN_HPP
