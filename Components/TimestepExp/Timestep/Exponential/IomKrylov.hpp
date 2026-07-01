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
#include "Timestep/Exponential/details/TimesteppperTools.hpp"

#include <iostream>
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
    * @param n    Size of distributed part of vector
    */
   int compute(Matrix& matV, Matrix& matH, const int j,
      const int m, const int n);

   /**
    * @brief Access action of A functor
    */
   TAfunc& aFunc();

   /**
    * @brief Tolerance
    */
   double tol() const;

   /**
    * @brief Gram-schmidt order
    */
   int p() const;

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
   const int jIn, const int m, const int n)
{
   // Check H is big enough
   assert(matH.rows() == matH.cols());
   assert(matH.rows() >= m+1);

   auto&& A = *this->mpAfunc;

   int j = jIn;
   for(; j < m; j++)
   {
      if(QuICCEnv().allowsIO())
      {
         std::cerr << "      - compute Jacobian" << std::endl;
      }

      // Build next vector
      A(matV.col(j+1), matV.col(j));

      // Orthogonalization
      int i0 = std::max(0, j + 1 - this->mcP);
      Matrix colH = details::computeAugmentedDot(matV, i0, j, matV, j+1, n);
      for(int i = i0; i <= j; i++)
      {
         matH(i, j) = colH(i-i0, 0);

         matV.col(j+1) -= matH(i,j)*matV.col(i);
      }

      // Norm
      double normV = details::computeAugmented2Norm(matV, j+1, n);
      if(QuICCEnv().allowsIO())
      {
         std::cerr << "         normV = " << normV << std::endl;
      }

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

template <typename TAfunc>
double IomKrylov<TAfunc>::tol() const
{
   return this->mcTol;
}

template <typename TAfunc>
int IomKrylov<TAfunc>::p() const
{
   return this->mcP;
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_IOMKRYLOV_HPP
