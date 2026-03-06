/**
 * @file ArnoldiKrylov.hpp
 * @brief Arnoldi iteration for building Krylov subspace
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_ARNOLDIKRYLOV_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_ARNOLDIKRYLOV_HPP

// System includes
//
#include <memory>

// Project includes
//

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

   /**
    * @brief Arnoldi iteration for building Krylov subspace
    */
   template <typename TAfunc>
   class ArnoldiKrylov
   {
      /**
       * @brief ctor
       */
      ArnoldiKrylov(std::unique_ptr<TAfunc>&& a);
      
      /**
       * @brief dtor
       */
      virtual ~ArnoldiKrylov() = default;

      /**
       * @brief Compute Krylov subspace approximation
       *
       * @param matV Krylov subspace basis
       * @param matH Projection of A on Krylov subspace
       * @param matB B matrix
       * @param j    starting index
       * @param m    Max size of Krylov subspace
       */
      int compute(matV, matH, const matB, const int j, const int m);

      private:
   /*
    * @brief Functor for action of A matrix
    */
   std::unique_ptr<TAfunc> mpAfunc;
   };

   template <typename TAfunc>
      ArnoldiKrylov<TAfunc>::ArnoldiKrylov(std::unique_ptr<TAfunc>&& a)
      : mpAfunc(std::move(a))
      {
      }

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_ARNOLDIKRYLOV_HPP
