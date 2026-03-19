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
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

/**
 * @brief Incomplete orthogonalization method for building Krylov subspace
 */
template <typename TAfunc> class IomKrylov
{
public:
   /**
    * @brief ctor
    *
    * @param a    Functor for action of A matrix
    * @param p    Length of incomplete orthogonalization
    * @param tol  Tolerance for subspace convergence
    */
   IomKrylov(std::shared_ptr<TAfunc>& a, const int p, const double tol);

   /**
    * @brief ctor
    *
    * @param a Functor for action of A matrix
    */
   IomKrylov(std::shared_ptr<TAfunc>& a);

   /**
    * @brief dtor
    */
   virtual ~IomKrylov() = default;

   /**
    * @brief Compute Krylov subspace approximation
    *
    * @param matV Krylov subspace basis
    * @param matH Projection of A on Krylov subspace
    * @param j    starting index
    * @param m    Max size of Krylov subspace
    */
   int compute(Matrix& matV, Matrix& matH, const int j,
      const int m);

   /**
    * @brief Access action of A functor
    */
   TAfunc& aFunc();

private:
   /**
    * @brief Length of incomplete orthogonalization
    */
   const int mcP;

   /**
    * @brief Convergence tolerance
    */
   const double mcTol;

   /*
    * @brief Functor for action of A matrix
    */
   std::shared_ptr<TAfunc> mpAfunc;
};

template <typename TAfunc>
IomKrylov<TAfunc>::IomKrylov(std::shared_ptr<TAfunc>& a, const int p, const double tol) :
    mcP(p), mcTol(tol), mpAfunc(a)
{}

template <typename TAfunc>
IomKrylov<TAfunc>::IomKrylov(std::shared_ptr<TAfunc>& a) : IomKrylov(std::move(a), 2, 1e-12)
{}

template <typename TAfunc>
int IomKrylov<TAfunc>::compute(Matrix& matV, Matrix& matH,
   const int jIn, const int m)
{
   // Check H is big enough
   assert(matH.rows() == matH.cols());
   assert(matH.rows() >= m+1);

   auto&& A = *this->mpAfunc;

   int j = jIn;
   for(; j < m; j++)
   {
      // Build next vector
      A(matV.col(j+1), matV.col(j));

      // Orthogonalization
      for(int i = std::max(0, j + 1 - this->mcP); i <= j; i++)
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
         j++;
         break;
      }

      matH(j+1, j) = normV;
      matV.col(j+1).array() /= normV;
   }

   return j;
}

template <typename TAfunc>
TAfunc& IomKrylov<TAfunc>::aFunc()
{
   return *this->mpAfunc;
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_IOMKRYLOV_HPP
