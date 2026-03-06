/**
 * @file IomKrylov.hpp
 * @brief Incomplete orthogonalization method for building Krylov subspace
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_IOMKRYLOV_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_IOMKRYLOV_HPP

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
    * @brief Incomplete orthogonalization method for building Krylov subspace
    */
   template <typename TAfunc>
   class IomKrylov
   {
      public:
      /**
       * @brief ctor
       *
       * @param a Functor for action of A matrix
       * @param p Length of incomplete orthogonalization
       */
      IomKrylov(std::unique_ptr<TAfunc>&& a, const int p);
      
      /**
       * @brief ctor
       *
       * @param a Functor for action of A matrix
       */
      IomKrylov(std::unique_ptr<TAfunc>&& a);
      
      /**
       * @brief dtor
       */
      virtual ~IomKrylov() = default;

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
      /**
       * @brief Length of incomplete orthogonalization
       */
      const int mcP;
   /*
    * @brief Functor for action of A matrix
    */
   std::unique_ptr<TAfunc> mpAfunc;
   };

   template <typename TAfunc>
      IomKrylov<TAfunc>::IomKrylov(std::unique_ptr<TAfunc>&& a, const int p)
      : mcP(p), mpAfunc(std::move(a))
      {
      }

   template <typename TAfunc>
      IomKrylov<TAfunc>::IomKrylov(std::unique_ptr<TAfunc>&& a)
      : IomKrylov(a, 2)
      {
      }

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_IOMKRYLOV_HPP
