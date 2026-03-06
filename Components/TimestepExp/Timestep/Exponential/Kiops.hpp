/**
 * @file Kiops.hpp
 * @brief KIOPS algorithm for exponential integrators
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP

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
 * @brief KIOPS algorithm for exponential integrators
 */
template <typename TKrylov, typename TExponential> class Kiops
{
public:
   /**
    * @brief ctor
    */
   Kiops(std::unique_ptr<TKrylov>&& K,
      std::unique_ptr<TExponential> E, const double tol, const int mMin, const int mMax);

   /**
    * @brief ctor
    */
   Kiops(std::unique_ptr<TKrylov>&& K,
      std::unique_ptr<TExponential> E);

   /**
    * @brief dtor
    */
   virtual ~Kiops() = default;

private:
   /**
    * @brief Tolerance
    */
   const double mcTol;

   /**
    * @brief Min size Krylov subspace
    */
   const int mcMMin;

   /**
    * @brief Max size of Krylov subspace
    */
   const int mcMMax;

   /*
    * @brief Algorithm for computing Krylov subspace
    */
   std::unique_ptr<TKrylov> mpKfunc;

   /*
    * @brief Algorithm for computing small dense matrix exponential
    */
   std::unique_ptr<TExponential> mpEfunc;

   /**
    * @brief Size of Krylov subspace
    */
   double mM;
};

template <typename TKrylov, typename TExponential>
Kiops<TKrylov, TExponential>::Kiops(std::unique_ptr<TKrylov>&& k, std::unique_ptr<TExponential> e, const double tol, const int mMin, const int mMax) :
    mcTol(tol), mcMMin(mMin), mcMMax(mMax), mpKfunc(std::move(k)), mpEfunc(std::move(e))
{}

template <typename TKrylov, typename TExponential>
Kiops<TKrylov, TExponential>::Kiops(std::unique_ptr<TKrylov>&& k, std::unique_ptr<TExponential> e) :
    Kiops(k, e, 1e-7, 10, 128)
{}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP
