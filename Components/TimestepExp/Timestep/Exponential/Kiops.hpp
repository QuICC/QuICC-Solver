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

// Project includes
//
#include "Types/Typedefs.hpp"

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

   /** 
    * @brief Compute action of phi
    */
   int compute(Matrix& matW, const std::vector<double>& ts, const Matrix& matU, const int m, const Task task) const;

private:
   struct Stats
   {
      int step = 0;
      int krystep = 0;
      int reject = 0;
      int exps = 0;
      int m_ret = 0;
      double conv = 0.0;
   };
   /**
    * @brief Initialize Krylov subspace
    */
   double initKrylov(Matrix& matV, const Matrix& matW, const double tNow, const int l, const int p) const;

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

template <typename TKrylov, typename TExponential>
int Kiops<TKrylov, TExponential>::compute(Matrix& matW, const std::vector<double>& ts, const Matrix& matU, const int mInit, const Task task) const
{
   const double sgn = std::copysign(1.0, ts.back());
   const double tEnd = std::abs(ts.back());
   const int nStep = ts.size();

   // Setting the safety factors and tolerance requirements
   double gammaMMax;
   double gamma;
   if(tEnd > 1)
   {
      gamma = 0.2;
      gammaMMax = 0.1;
   }
   else
   {
      gamma = 0.9;
      gammaMMax = 0.6;
   }
   const double delta = 1.4;

   // Get dimensions
   const int n = matU.rows();
   const int p = matU.cols() - 1;

   double tNow = 0;
   double tau = tEnd;

   bool happy = false;

   // Limit Krylov subspace size
   int m = std::max(this->mcMMin, std::min(mInit, this->mcMMax));

   int j = 0;
   int l = 0;
   int ireject = 0;

   Stats stats;

   // Initialize variables
   matW.resize(n, nStep);
   Matrix matH = Matrix::Zero(m + 1, m + 1);
   Matrix matV = Matrix::Zero(n + p, m + 1);
   
   // Adaptive part
   int oldm = -1;
   double oldTau = 0.0;
   double omega = 0.0;
   bool orderOld = true;
   bool kestOld = true;

   double order = 0.0;
   double beta = 0.0;
   double kest = 0.0;

   matW.col(0) = matU.col(0);
   while( tNow < ts.back() )
   {
      // Initial vector for Krylov
      if(j == 0)
      {
         beta = this->initKrylov(matV, matW, tNow, l, p);

         matH.setZero();
      }

      // Compute Krylov subspace
      j = this->mpKfunc->compute(matV, matH, j, m);
      happy = (j <= m);
      matH(0, j) = 1.0;
      double hm1m = matH(j,j-1);
      matH(j, j-1) = 0.0;

      // Compute exponential of H
      Matrix matF = sgn * tau * matH.block(0, 0, j + 1, j + 1);
      this->mpEfunc(matF);
      stats.exps++;

      // Restore H
      matH(j, j-1) = hm1m;

      double err;
      double tNew;
      int mNew;
      double tauNew;

      // Krylov subspace converged before m
      if(happy)
      {
         omega = 0.0;
         err = 0.0;
         tNew = std::min(tEnd - (tNow + tau), tau);
         happy = false;
      }
      else
      {
         // Local truncation error
         err = std::abs(beta * hm1m * matF(j - 1, j));

         // Error for this step
         double oldOmega = omega;
         omega = tEnd * err / (tau * this->mcTol);

         // Estimate order
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

         // Estimate k
         if(m != oldm && tau == oldTau && ireject >= 1)
         {
            kest = std::max(1.1, std::pow((omega / oldOmega), 1./(oldm - m)));
            kestOld = true;
         }
         else if(kestOld || ireject == 0)
         {
            kestOld = false;
            kest = 2.0;
         }
         else
         {
            kestOld = true;
         }

         double remainingTime;
         if(omega > delta)
         {
            remainingTime = tEnd - tNow;
         }
         else
         {
            remainingTime = tEnd - (tNow + tau);
         }
         
         // Adaptive Krylov
         double sameTau = std::min(remainingTime, tau);
         double tOpt = tau * std::pow((gamma / omega), 1./order);
         tOpt = std::min(remainingTime, std::max(tau / 5., std::min(5. * tau, tOpt)));

         int mOpt = std::ceil(j + std::log(omega / gamma) / std::log(kest));
         mOpt = std::max(this->mcMMin, std::min(this->mcMMax, std::max(static_cast<int>(std::floor(3./4. * m)), std::min(mOpt, static_cast<int>(std::ceil(4./3. * m))))));

         if( j == this->mcMMax)
         {
            if(omega > delta)
            {
               mNew = j;
               tNew = tau * std::pow(gammaMMax / omega, 1./order);
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

      // Achieved required tolerance
      if( omega <= delta)
      {
         stats.reject += ireject;
         stats.step += 1;

         int blownTs = 0;
         double nextT = tNow + tau;
         for(int k = l; k < nStep; k++)
         {
            if(std::abs(ts.at(k)) < std::abs(nextT))
            {
               blownTs += 1;
            }
         }

         if(blownTs != 0)
         {
            matW.col(l + blownTs) = matW.col(l);

            for(int k = 0; k < blownTs; k++)
            {
               double tPhantom = ts.at(l + k) - tNow;
               Matrix matF2 = sgn * tPhantom * matH.block(0, 0, j, j);
               this->mpEfunc(matF2);
               matW.col(l + k) = beta * matV.block(0, 0, n, j) * matF2.col(0);
            }
            l += blownTs;
         }

         Matrix tmp1 = beta * matF.col(0);
         Matrix tmp2 = matW.col(l);
         tmp2 = tmp1(0) * matV.block(0,0,n,1);
         for(int i = 1; i < j; i++)
         {
            tmp2 = tmp1(i) * matV.block(0,i,n,1);
         }
         tNow += tau;

         j = 0;
         ireject = 0;

         stats.conv += err;
      }
      else
      {
         ireject += 1;

         // Restore H
         matH(0, j) = 0.0;
      }

      oldTau = tau;
      tau = tauNew;

      oldm = m;
      m = mNew;
   }

   if(task == Task::I)
   {
      for(int k = 0; k < nStep; k++)
      {
         matW.col(k).array() /= std::pow(ts.at(k), p);
      }
   }

   stats.m_ret = m;

   return m;
}

template <typename TKrylov, typename TExponential>
double Kiops<TKrylov, TExponential>::initKrylov(Matrix& matV, const Matrix& matW, const double tNow, const int l, const int p) const
{
   matV.col(0) = matW.col(l);

   int iLast = matV.rows() - 1;
   matV(iLast, l) = 1;
   int pf = 1;
   for(int i = 1; i < p; i++)
   {
      pf *= i;
      matV(iLast - i, l) = std::pow(tNow, i)/static_cast<double>(pf);
   }

   double beta = matV.col(0).squaredNorm();
   // MPI HERE
   beta = std::sqrt(beta);
   matV.col(0).array() /= beta;

   return beta;
}

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP
