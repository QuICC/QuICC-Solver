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
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

/**
 * @brief Arnoldi iteration for building Krylov subspace
 */
template <typename TAfunc> class ArnoldiKrylov
{
public:
   /**
    * @brief ctor
    *
    * @param a    Functor for action of A matrix
    * @param tol  Tolerance for subspace convergence
    */
   ArnoldiKrylov(std::unique_ptr<TAfunc>&& a, const double tol);

   /**
    * @brief ctor
    *
    * @param a    Functor for action of A matrix
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
    * @param j    starting index
    * @param m    Max size of Krylov subspace
    */
   int compute(Matrix& matV, Matrix& matH, const int j, const int m);

   /**
    * @brief Access action of A functor
    */
   TAfunc& aFunc();

private:
   /**
    * @brief Convergence tolerance
    */
   const double mcTol;

   /*
    * @brief Functor for action of A matrix
    */
   std::unique_ptr<TAfunc> mpAfunc;
}
;

template <typename TAfunc>
ArnoldiKrylov<TAfunc>::ArnoldiKrylov(std::unique_ptr<TAfunc>&& a,
   const double tol) :
    mcTol(tol), mpAfunc(std::move(a))
{}

template <typename TAfunc>
ArnoldiKrylov<TAfunc>::ArnoldiKrylov(std::unique_ptr<TAfunc>&& a) :
    ArnoldiKrylov(std::move(a), 1e-12)
{}

template <typename TAfunc>
int ArnoldiKrylov<TAfunc>::compute(Matrix& matV, Matrix& matH, const int jIn,
   const int m)
{
   // Check H is big enough
   assert(matH.rows() == matH.cols());
   assert(matH.rows() >= m+1);
   assert(jIn < m);

   auto&& A = *this->mpAfunc;

   int j = jIn;
   for(; j < m; j++)
   {
      // Build next vector
      A(matV.col(j+1), matV.col(j));

      // Orthogonalization
      for(int i = 0; i <= j; i++)
      {
         matH(i, j) = matV.col(i).dot(matV.col(j+1));
         // MPI VERSION HERE
         matV.col(j+1) -= matH(i,j)*matV.col(i);
      }

      // Norm
      double normV = matV.col(j+1).squaredNorm();
      // MPI VERSION HERE
      normV = std::sqrt(normV);

      // Stop if subspace converged sufficiently
      if(normV < this->mcTol)
      {
         break;
      }

      matH(j+1, j) = normV;
      matV.col(j+1).array() /= normV;
   }

   return j;
}

template <typename TAfunc>
TAfunc& ArnoldiKrylov<TAfunc>::aFunc()
{
   return *this->mpAfunc;
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_ARNOLDIKRYLOV_HPP
