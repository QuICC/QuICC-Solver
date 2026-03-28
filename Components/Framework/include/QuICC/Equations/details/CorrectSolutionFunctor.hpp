/**
 * @file CorrectSolutionFunctor.hpp
 * @brief Implementation of the CorrectSolution functors
 */

#ifndef QUICC_EQUATIONS_DETAILS_CORRECTSOLUTIONFUNCTOR_HPP
#define QUICC_EQUATIONS_DETAILS_CORRECTSOLUTIONFUNCTOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Equations {

namespace details {

template <CouplingIndexType IndexType> class CorrectSolutionFunctor
{
public:
   /**
    * @brief ctor
    */
   CorrectSolutionFunctor(const Resolution& res, const CouplingInformation& cinfo, const std::vector<std::tuple<MHDVariant, int, int, int>>& corrections,
      const int matIdx);

   /**
    * @brief deleted default ctor
    */
   CorrectSolutionFunctor() = delete;

   /**
    * @brief dtor
    */
   ~CorrectSolutionFunctor() = default;

   /**
    * @brief Apply functor
    */
   template <typename TData> void apply(TData& storage, const int start);

   /**
    * Brief Corrections
    */
   std::vector<std::tuple<MHDComplex, int, int>> corr;

private:
   /**
    * @brief Init functor
    */
   void init(
      const std::vector<std::tuple<MHDVariant, int, int, int>>& corrections);

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
};

template <CouplingIndexType IndexType>
CorrectSolutionFunctor<IndexType>::CorrectSolutionFunctor(
   const Resolution& res, const CouplingInformation& cinfo,
   const std::vector<std::tuple<MHDVariant, int, int, int>>& corrections,
   const int matIdx) :
    res(res), cinfo(cinfo), matIdx(matIdx)
{
   this->init(corrections);
}

} // namespace details
} // namespace Equations
} // namespace QuICC

#endif // QUICC_EQUATIONS_DETAILS_CORRECTSOLUTIONFUNCTOR_HPP
