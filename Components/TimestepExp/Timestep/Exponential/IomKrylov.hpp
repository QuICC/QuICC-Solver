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
template <typename TAfunc, typename TOfunc> class IomKrylov
{
public:
   /**
    * @brief ctor
    *
    * @param a    Functor for action of A matrix
    * @param o    Functor for orthogonalization
    * @param tol  Tolerance for subspace convergence
    */
   IomKrylov(std::shared_ptr<TAfunc>& a, std::shared_ptr<TOfunc>& o, const double tol);

   /**
    * @brief ctor
    *
    * @param a Functor for action of A matrix
    */
   IomKrylov(std::shared_ptr<TAfunc>& a, std::shared_ptr<TOfunc>& o);

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
    * @brief Access action of A functor
    */
   TOfunc& oFunc();

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
    * @brief Convergence tolerance
    */
   const double mcTol;

   /*
    * @brief Functor for action of A matrix
    */
   std::shared_ptr<TAfunc> mpAfunc;

   /**
    * @brief Functor for orthogonalization
    */
   std::shared_ptr<TOfunc> mpOfunc;
};

template <typename TAfunc, typename TOfunc>
IomKrylov<TAfunc,TOfunc>::IomKrylov(std::shared_ptr<TAfunc>& a, std::shared_ptr<TOfunc>& o, const double tol) :
    mcTol(tol), mpAfunc(a), mpOfunc(o)
{}

template <typename TAfunc, typename TOfunc>
IomKrylov<TAfunc,TOfunc>::IomKrylov(std::shared_ptr<TAfunc>& a, std::shared_ptr<TOfunc>& o) : IomKrylov(std::move(a), std::move(o), 1e-12)
{}

template <typename TAfunc, typename TOfunc>
int IomKrylov<TAfunc,TOfunc>::compute(Matrix& matV, Matrix& matH,
   const int jIn, const int m, const int n)
{
   // Check H is big enough
   assert(matH.rows() == matH.cols());
   assert(matH.rows() >= m+1);

   auto&& A = *this->mpAfunc;
   auto&& ortho = *this->mpOfunc;

   int j = jIn;
   for(; j < m; j++)
   {
      // Build next vector
      A(matV.col(j+1), matV.col(j));

      // Orthogonalize
      double normV = ortho(matV, matH, j, n);

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

template <typename TAfunc, typename TOfunc>
TAfunc& IomKrylov<TAfunc,TOfunc>::aFunc()
{
   return *this->mpAfunc;
}

template <typename TAfunc, typename TOfunc>
TOfunc& IomKrylov<TAfunc,TOfunc>::oFunc()
{
   return *this->mpOfunc;
}

template <typename TAfunc, typename TOfunc>
double IomKrylov<TAfunc,TOfunc>::tol() const
{
   return this->mcTol;
}

template <typename TAfunc, typename TOfunc>
int IomKrylov<TAfunc,TOfunc>::p() const
{
   return this->mpOfunc->p();
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_IOMKRYLOV_HPP
