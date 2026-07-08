/**
 * @file Kiops.hpp
 * @brief KIOPS algorithm for exponential integrators
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP

// System includes
//
#include <cmath>
#include <memory>
#include <iostream>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "Types/Typedefs.hpp"
#include "Timestep/Exponential/details/TimesteppperTools.hpp"
#include "Timestep/Exponential/RuntimeStatistics.hpp"

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
   enum class Task {I, II};

   /**
    * @brief ctor
    */
   Kiops(std::unique_ptr<TKrylov>&& K,
      std::unique_ptr<TExponential> E, const double tol, const double delta, const int mMin, const int mMax);

   /**
    * @brief ctor
    */
   Kiops(std::unique_ptr<TKrylov>&& K,
      std::unique_ptr<TExponential> E);

   /**
    * @brief dtor
    */
   virtual ~Kiops() = default;

   /**
    * @brief setup
    */
   void setup(const std::vector<double>& ts, const Matrix& matU);

   /**
    * @brief Compute action of phi
    */
   int compute(Matrix& matW, const std::vector<double>& ts, const Matrix& matU, const int m, const Task task) const;

   /**
    * @brief Access Krylov subspace functor
    */
   TKrylov& kFunc();

   /**
    * @brief Access explicit dense exponential functor
    */
   TExponential& eFunc();

   /**
    * @brief Print runtime information
    */
   void printInfo() const;

private:
   /**
    * @brief Scaling factors for B
    */
   std::pair<double,double> computeNuMu(const MHDFloat norm) const;

   /**
    * @brief Initialize Krylov subspace
    */
   double initKrylov(Matrix& matV, const Matrix& matW, const double tNow, const int l, const int p, const std::pair<double,double>& numu) const;

   /**
    * @brief Estimate order
    */
   void estimateOrder(double &order, bool &orderOld, const double omega, const double omegaOld, const double tau, const double oldTau, const int m, const int oldm, const int j, const int ireject) const;

   /**
    * @brief Estimate k
    */
   void estimateK(double &kest, bool& kestOld, const double omega, const double omegaOld, const double tau, const double oldTau, const int m, const int oldm, const int ireject) const;

   /**
    * @brief Adaptive Krylov subspace
    */
   void adaptiveKrylov(int& mNew, double& tNew, double& omega, bool& orderOld, bool& kestOld, const double err, const double tNow, const double tEnd, const double tau, const double oldTau, const int m, const int oldm,  const int j, const int ireject) const;

   /**
    * @brief Update solution
    */
   void updateSolution(int& l, double& tNow, Matrix& matW, const double tau, const double beta, const std::vector<double>& ts, const Matrix& matV, const Matrix& matH, const Matrix& matF, const int n, const int j) const;

   /**
    * @brief Define action of augmented A
    */
   void defineAugmentedA(const Matrix& matU);

   /**
    * @brief Set safety factors
    */
   void setSafety(const std::vector<double>& ts);

   /**
    * @brief Tolerance
    */
   const double mcTol;

   /**
    * @brief Scaled error
    */
   const double mcDelta;

   /**
    * @brief Min size Krylov subspace
    */
   const int mcMMin;

   /**
    * @brief Max size of Krylov subspace
    */
   const int mcMMax;

   /**
    * @brief Safety factor for scaled error (first), if Mmax is reached(second)
    */
   std::pair<double,double> mGamma;

   /**
    * @brief Scaling factors Nu and Mu from U
    */
   std::pair<double, double> mNuMu;

   /*
    * @brief Algorithm for computing Krylov subspace
    */
   std::unique_ptr<TKrylov> mpKfunc;

   /*
    * @brief Algorithm for computing small dense matrix exponential
    */
   std::unique_ptr<TExponential> mpEfunc;

   /**
    * @brief Runtime statistatics
    */
   mutable RuntimeStatistics mStats;
};

template <typename TKrylov, typename TExponential>
Kiops<TKrylov, TExponential>::Kiops(std::unique_ptr<TKrylov>&& k, std::unique_ptr<TExponential> e, const double tol, const double delta, const int mMin, const int mMax) :
    mcTol(tol), mcDelta(delta), mcMMin(mMin), mcMMax(mMax), mpKfunc(std::move(k)), mpEfunc(std::move(e)), mStats()
{}

template <typename TKrylov, typename TExponential>
Kiops<TKrylov, TExponential>::Kiops(std::unique_ptr<TKrylov>&& k, std::unique_ptr<TExponential> e) :
    Kiops(std::move(k), std::move(e), 1e-12, 1.4, 10, 128)
{}

template <typename TKrylov, typename TExponential>
void Kiops<TKrylov, TExponential>::setup(const std::vector<double>& ts, const Matrix& matU)
{
   this->setSafety(ts);

   this->defineAugmentedA(matU);
}

template <typename TKrylov, typename TExponential>
void Kiops<TKrylov, TExponential>::defineAugmentedA(const Matrix& matU)
{
   Matrix matB;
   MHDFloat normB;
   if(matU.cols() > 1)
   {
      matB = matU.rightCols(matU.cols()-1).rowwise().reverse();
      normB = details::compute1Norm(matB);
   }
   else
   {
      matB = Matrix::Zero(matU.rows(), 1);
      normB = 0;
   }
   this->mNuMu = this->computeNuMu(normB);

   if(matU.cols() > 1)
   {
      const double& nu = this->mNuMu.first;
      matB.array() *= nu;
   }

   this->mpKfunc->aFunc().updateB(matB);
}

template <typename TKrylov, typename TExponential>
int Kiops<TKrylov, TExponential>::compute(Matrix& matW, const std::vector<double>& ts, const Matrix& matU, const int mInit, const Task task) const
{
   // reset stats
   this->mStats.reset();
   this->mStats.gram_p = this->mpKfunc->p();

   const double sgn = std::copysign(1.0, ts.back());
   const double tEnd = std::abs(ts.back());

   // Get dimensions
   const int n = matU.rows();
   const int p = std::max(1, static_cast<int>(matU.cols()) - 1);

   double tNow = 0;
   double tau = tEnd;

   bool happy = false;

   // Limit Krylov subspace size
   int m = std::max(this->mcMMin, std::min(mInit, this->mcMMax));

   int j = 0;
   int l = 0;
   int ireject = 0;

   // Initialize variables
   matW.resize(n, ts.size());
   Matrix matH = Matrix::Zero(this->mcMMax + 1, this->mcMMax + 1);
   Matrix matV = Matrix::Zero(n + p, this->mcMMax + 1);

   // Adaptive part
   int oldm = -1;
   bool kestOld = true;
   bool orderOld = true;
   double oldTau = 0.0;
   double omega = 0.0;

   double beta = 0.0;

   matW.col(0) = matU.col(0);
   while( tNow < ts.back() )
   {
      if(QuICCEnv().allowsIO())
      {
         std::cerr << "   --- KIOPS iteration j = " << j << ", m = " << m << ", tau = " << tau << std::endl;
      }

      // Initial vector for Krylov
      if(j == 0)
      {
         beta = this->initKrylov(matV, matW, tNow, l, p, this->mNuMu);

         matH.setZero();
      }

      // Compute Krylov subspace
      j = this->mpKfunc->compute(matV, matH, j, m, n);
      happy = (j < m);
      matH(0, j) = 1.0;
      double hm1m = matH(j,j-1);
      matH(j, j-1) = 0.0;

      // Compute exponential of H
      Matrix matF = sgn * tau * matH.block(0, 0, j + 1, j + 1);
      matF = this->mpEfunc->compute(matF);
      this->mStats.exps++;

      // Restore H
      matH(j, j-1) = hm1m;

      double err;
      double tNew;
      int mNew;

      // Krylov subspace converged before m
      if(happy)
      {
         omega = 0.0;
         err = 0.0;
         tNew = std::min(tEnd - (tNow + tau), tau);
         mNew = j;
         happy = false;
      }
      else
      {
         // Local truncation error
         err = std::abs(beta * hm1m * matF(j - 1, j));

         // Adaptive Krylov
         this->adaptiveKrylov(mNew, tNew, omega, orderOld, kestOld, err, tNow, tEnd, tau, oldTau, m, oldm, j, ireject);
      }

      // Achieved required tolerance
      if(omega <= this->mcDelta)
      {
         this->mStats.reject += ireject;
         this->mStats.step += 1;

         this->updateSolution(l, tNow, matW, tau, beta, ts, matV, matH, matF, n, j);

         j = 0;
         ireject = 0;

         this->mStats.conv += err;
      }
      // Did not yet converge to tolerance
      else
      {
         ireject += 1;

         // Restore H
         matH(0, j) = 0.0;
      }

      oldTau = tau;
      tau = tNew;

      oldm = m;
      m = mNew;
   }

   // Rescale by 1/t^p for Task I
   if(task == Task::I)
   {
      for(std::size_t k = 0; k < ts.size(); k++)
      {
         matW.col(k).array() /= std::pow(ts.at(k), p);
      }
   }

   this->mStats.m = m;
   this->mStats.update();

   return m;
}

template <typename TKrylov, typename TExponential>
double Kiops<TKrylov, TExponential>::initKrylov(Matrix& matV, const Matrix& matW, const double tNow, const int l, const int p, const std::pair<double,double>& numu) const
{
   const double& mu = numu.second;

   matV.col(0).topRows(matW.rows()) = matW.col(l);

   int iLast = matV.rows() - 1;
   matV(iLast, l) = mu;
   int pf = 1;
   for(int i = 1; i < p; i++)
   {
      pf *= i;
      matV(iLast - i, 0) = std::pow(tNow, i)*mu/static_cast<double>(pf);
   }

   double beta = details::computeAugmented2Norm(matV, 0, matW.rows());
   matV.col(0).array() /= beta;

   return beta;
}

template <typename TKrylov, typename TExponential>
void Kiops<TKrylov, TExponential>::estimateOrder(double& order, bool& orderOld, const double omega, const double oldOmega, const double tau, const double oldTau, const int m, const int oldm, const int j, const int ireject) const
{
   if(m == oldm && tau != oldTau && ireject >= 1)
   {
      order = std::max(1.0, std::log(omega/oldOmega) / std::log(tau/oldTau));
      orderOld = false;
   }
   else if(orderOld || ireject == 0)
   {
      orderOld = true;
      order = j / 4.;
   }
   else
   {
      orderOld = true;
   }
}

template <typename TKrylov, typename TExponential>
void Kiops<TKrylov, TExponential>::estimateK(double &kest, bool& kestOld, const double omega, const double oldOmega, const double tau, const double oldTau, const int m, const int oldm, const int ireject) const
{
   if(m != oldm && tau == oldTau && ireject >= 1)
   {
      kest = std::max(1.1, std::pow((omega / oldOmega), 1./(oldm - m)));
      kestOld = false;
   }
   else if(kestOld || ireject == 0)
   {
      kestOld = true;
      kest = 2.0;
   }
   else
   {
      kestOld = true;
   }
}

template <typename TKrylov, typename TExponential>
void Kiops<TKrylov, TExponential>::adaptiveKrylov(int& mNew, double& tNew, double& omega, bool& orderOld, bool& kestOld, const double err, const double tNow, const double tEnd, const double tau, const double oldTau, const int m, const int oldm,  const int j, const int ireject) const
{
   double order = 0.0;
   double kest = 0.0;

   // Error for this step
   double oldOmega = omega;
   omega = tEnd * err / (tau * this->mcTol);

   if(std::isnan(omega))
   {
      throw std::logic_error("Error estimate is NaN");
   }

   // Estimate order: order + 1 = 1/(q + 1) in paper
   this->estimateOrder(order, orderOld, omega, oldOmega, tau, oldTau, m, oldm, j, ireject);
   // Estimate k
   this->estimateK(kest, kestOld, omega, oldOmega, tau, oldTau, m, oldm, ireject);

   double remainingTime;
   if(omega > this->mcDelta)
   {
      remainingTime = tEnd - tNow;
   }
   else
   {
      remainingTime = tEnd - (tNow + tau);
   }

   double sameTau = std::min(remainingTime, tau);
   double tOpt = tau * std::pow((this->mGamma.first / omega), order + 1);
   tOpt = std::min(remainingTime, std::max(tau / 5., std::min(5. * tau, tOpt)));

   int mOpt = std::ceil(j + std::log(omega / this->mGamma.first) / std::log(kest));
   if(QuICCEnv().allowsIO())
   {
      std::cerr << "      - Adaptive Krylovo: omega = " << omega << ", order = " << order << ", kest = " << kest << ", mOpt = " << mOpt << ", m = " << m  << std::endl;
   }
   mOpt = std::max(this->mcMMin, std::min(this->mcMMax, std::max(static_cast<int>(std::floor(3./4. * m)), std::min(mOpt, static_cast<int>(std::ceil(4./3. * m))))));

   if(j == this->mcMMax)
   {
      if(omega > this->mcDelta)
      {
         mNew = j;
         tNew = tau * std::pow(this->mGamma.second / omega, order + 1);
         tNew = std::min(tEnd - tNow, std::max(tau / 5., tNew));
      }
      else
      {
         tNew = tOpt;
         mNew = m;
      }
   }
   else
   {
      mNew = mOpt;
      tNew = sameTau;
   }
}

template <typename TKrylov, typename TExponential>
void Kiops<TKrylov, TExponential>::updateSolution(int& l, double& tNow, Matrix& matW, const double tau, const double beta, const std::vector<double>& ts, const Matrix& matV, const Matrix& matH, const Matrix& matF, const int n, const int j) const
{
   const double sgn = std::copysign(1.0, ts.back());

   int nTau = 0;
   double tNext = tNow + tau;
   for(std::size_t k = l; k < ts.size(); k++)
   {
      if(std::abs(ts.at(k)) < std::abs(tNext))
      {
         nTau++;
      }
   }

   if(nTau != 0)
   {
      matW.col(l + nTau) = matW.col(l);

      for(int k = 0; k < nTau; k++)
      {
         double tbar = ts.at(l + k) - tNow;
         Matrix matF2 = sgn * tbar * matH.block(0, 0, j, j);
         matF2 = this->mpEfunc->compute(matF2);
         matW.col(l + k) = beta * matV.block(0, 0, n, j) * matF2.col(0);
      }
      l += nTau;
   }

   Matrix tmp1 = beta * matF.col(0);
   matW.col(l) = tmp1(0) * matV.block(0,0,n,1);
   for(int i = 1; i < j; i++)
   {
      matW.col(l) += tmp1(i) * matV.block(0,i,n,1);
   }
   tNow += tau;
}

template <typename TKrylov, typename TExponential>
void Kiops<TKrylov, TExponential>::setSafety(const std::vector<double>& ts)
{
   const double tEnd = std::abs(ts.back());

   // Setting the safety factors and tolerance requirements
   if(tEnd > 1)
   {
      this->mGamma.first = 0.2;
      this->mGamma.second = 0.1;
   }
   else
   {
      this->mGamma.first = 0.9;
      this->mGamma.second = 0.6;
   }
}

template <typename TKrylov, typename TExponential>
std::pair<double,double> Kiops<TKrylov, TExponential>::computeNuMu(const MHDFloat normB) const
{
   double nu, mu;

   if(normB > 0)
   {
      double ex = std::ceil(std::log2(normB));
      nu = std::pow(2.0, -ex);
      mu = std::pow(2.0, ex);
   }
   else
   {
      nu = 1.0;
      mu = 1.0;
   }

   return std::make_pair(nu, mu);
}

template <typename TKrylov, typename TExponential>
TKrylov& Kiops<TKrylov, TExponential>::kFunc()
{
   return *this->mpKfunc;
}

template <typename TKrylov, typename TExponential>
TExponential& Kiops<TKrylov, TExponential>::eFunc()
{
   return *this->mpEfunc;
}

template <typename TKrylov, typename TExponential>
void Kiops<TKrylov, TExponential>::printInfo() const
{
   if(QuICCEnv().allowsIO())
   {
      std::cerr
         << std::string(5, '=')
         << " KIOPS information "
         << std::string(5, '=')
         << std::endl;

      this->mStats.printInfo();

      std::cerr
         << std::string(29, '=')
         << std::endl;
   }
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP
