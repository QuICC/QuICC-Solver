/**
 * @file TransposeFunctor.hpp
 * @brief Functor to compute data transpose between spatial schemes
 */

#ifndef QUICC_SOLVER_DETAILS_TRANSPOSEFUNCTOR_HPP
#define QUICC_SOLVER_DETAILS_TRANSPOSEFUNCTOR_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "ViewOps/Transpose/OpGrouped.hpp"
#include "Environment/QuICCEnv.hpp"

namespace QuICC {

namespace Solver {

namespace details {

template <typename T> class TransposeFunctor
{
public:
   /**
    * @brief ctor
    */
   TransposeFunctor() = default;

   /**
    * @brief dtor
    */
   ~TransposeFunctor() = default;

   /**
    * @brief Compute and add the explicit linear terms
    *
    * @param rSolverField  Solver field values
    * @param eqStart       Start index for the equation field
    * @param explicitField Explicit linear field values
    */
   const Framework::Selector::ScalarField<T>* apply(const Resolution& eqRes, const Framework::Selector::ScalarField<T>& field);

private:
   using VoutTy = View::View<T, View::DCCSC3D>;
   using VinTy = View::View<T, View::DCCSC3D>;
#ifdef QUICC_MPI
   using TransposeOpType = QuICC::Transpose::Mpi::OpGrouped<std::vector<VoutTy>, std::vector<VinTy>, Transpose::p021_t>;
#else
   using TransposeOpType = QuICC::Transpose::Cpu::OpGrouped<std::vector<VoutTy>, std::vector<VinTy>, Transpose::p021_t>;
#endif

   /**
    * @brief Temporary real storage
    */
   std::shared_ptr<Framework::Selector::ScalarField<T>> mspTmp;

   /**
    * @brief Transpose operator
    */
   std::unique_ptr<TransposeOpType> mTransOp;
};

template <typename T> const Framework::Selector::ScalarField<T>* TransposeFunctor<T>::apply(const Resolution& eqRes, const Framework::Selector::ScalarField<T>& field)
{
   QuICCEnv().synchronize();

   // Initialize data and operator
   if(!this->mspTmp)
   {
      auto spSetup = eqRes.spSpectralSetup();
      this->mspTmp = std::make_shared<Framework::Selector::ScalarField<T>>(spSetup);

      // Create transpose operator
#ifdef QUICC_MPI
      this->mTransOp = std::make_unique<TransposeOpType>(spSetup->mem());
#else
      this->mTransOp = std::make_unique<TransposeOpType>();
#endif
   }

   // Apply transpose
   std::vector<VoutTy> viewsOut = {this->mspTmp->rGlobalView()};
   std::vector<VinTy> viewsIn = {field.globalView()};
   this->mTransOp->apply(viewsOut, viewsIn);

   QuICCEnv().synchronize();

   return this->mspTmp.get();
}

} // namespace details
} // namespace Solver
} // namespace QuICC

#endif // QUICC_SOLVER_DETAILS_TRANSPOSEFUNCTOR_HPP
