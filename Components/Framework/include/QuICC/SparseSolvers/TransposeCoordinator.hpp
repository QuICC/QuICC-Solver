/**
 * @file TransposeCoordinator.hpp
 * @brief Coordinator for transpose operators
 */

#ifndef QUICC_SOLVER_TRANSPOSECOORDINATOR_HPP
#define QUICC_SOLVER_TRANSPOSECOORDINATOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SparseSolvers/details/TransposeFunctor.hpp"

namespace QuICC {

namespace Solver {

class TransposeCoordinator
{
public:
   /**
    * @brief ctor
    */
   TransposeCoordinator() = default;

   /**
    * @brief dtor
    */
   ~TransposeCoordinator() = default;

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param rSolverField  Solver field values
    * @param eqStart       Start index for the equation field
    * @param explicitField Explicit linear field values
    */
   template <typename T>
   const Framework::Selector::ScalarField<T>* apply(const Resolution& fieldRes, const Resolution& eqRes, const Framework::Selector::ScalarField<T>& field);

private:
   /**
    * @brief Get correct transpose operators
    */
   template <typename T>
      std::map<std::size_t, std::unique_ptr<details::TransposeFunctor<T>>>& transOp();

   /**
    * @brief Transpose functor for real field
    */
   std::map<std::size_t, std::unique_ptr<details::TransposeFunctor<MHDFloat>>> mRealTrans;

   /**
    * @brief Transpose functor for complex field
    */
   std::map<std::size_t, std::unique_ptr<details::TransposeFunctor<MHDComplex>>> mCplxTrans;
};

template <typename T>
   const Framework::Selector::ScalarField<T>* TransposeCoordinator::apply(const Resolution& fieldRes, const Resolution& eqRes, const Framework::Selector::ScalarField<T>& field)
{
   auto id = eqRes.sim().ss().id();
   if(id != fieldRes.sim().ss().id())
   {
      if(this->transOp<T>().count(id) == 0)
      {
         this->transOp<T>().emplace(id, std::make_unique<details::TransposeFunctor<T>>());
      }

      return this->transOp<T>().at(id)->apply(eqRes, field);
   }
   {
      return &field;
   }
}

template <typename T>
   std::map<std::size_t, std::unique_ptr<details::TransposeFunctor<T>>>& TransposeCoordinator::transOp()
{
   if constexpr(std::is_same_v<T, MHDFloat>)
   {
      return this->mRealTrans;
   }
   else if constexpr(std::is_same_v<T, MHDComplex>)
   {
      return this->mCplxTrans;
   }
   else
   {
      throw std::logic_error("Unknown field data for transpose");
   }
}

} // namespace Solver
} // namespace QuICC

#endif // QUICC_SOLVER_TRANSPOSECOORDINATOR_HPP
