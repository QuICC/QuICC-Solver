/**
 * @file ImExRK3.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order ~3
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Timestep/ImExRK3.hpp"

// Project includes
//

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   int ImExRK3::steps() const
   {
      return 3;
   }

   // Scheme order
   int ImExRK3::order() const
   {
      return 2;
   }

   // Scheme requires field at t_n
   int ImExRK3::fieldMemory() const
   {
      return 0;
   }

   // Scheme requires nonlinear term at t(n-1)
   int ImExRK3::nonlinearMemory() const
   {
      return 1;
   }

   // Name of the scheme
   std::string ImExRK3::name() const
   {
      return "ImExRK3";
   }

   ImExRK3::ImExRK3()
      : IImExOldScheme(),
      mAlpha(Eigen::Array<MHDFloat,3,1>(29./96., -3./40., 1./6.)),
      mBeta(Eigen::Array<MHDFloat,3,1>(37./160., 5./24., 1./6.)),
      mGamma(Eigen::Array<MHDFloat,3,1>(8./15., 5./12., 3./4.)),
      mZeta(Eigen::Array<MHDFloat,3,1>(0., -17./60., -5./12.)),
      mCEx(Eigen::Array<MHDFloat,3,1>(0., 0., 1.))
   {
   }

   ImExRK3::~ImExRK3()
   {
   }

   MHDFloat ImExRK3::lhsT(const int step) const
   {
      assert(step < this->steps());

      return 1./this->mBeta(step);
   }

   MHDFloat ImExRK3::lhsL(const int step) const
   {
      assert(step < this->steps());

      return 1.0;
   }

   MHDFloat ImExRK3::rhsT(const int i, const int step) const
   {
      assert(step < this->steps());
      assert(i > -1);
      assert(i < this->fieldMemory()+1);

      return 1./this->mBeta(step);
   }

   MHDFloat ImExRK3::rhsL(const int i, const int step) const
   {
      assert(step < this->steps());
      assert(i > -1);
      assert(i < this->fieldMemory()+1);

      return this->mAlpha(step)/this->mBeta(step);
   }

   MHDFloat ImExRK3::rhsN(const int i, const int step) const
   {
      assert(step < this->steps());
      assert(i > -1);
      assert(i < this->nonlinearMemory()+1);

      MHDFloat coeff;
      if(i == 0)
      {
         coeff = this->mGamma(step)/this->mBeta(step);
      } else
      {
         coeff = this->mZeta(step)/this->mBeta(step);
      }

      return coeff;
   }

   MHDFloat ImExRK3::cEx(const int i) const
   {
      return this->mCEx(i);
   }

   int ImExRK3::fieldMemory(const int step) const
   {
      assert(step < this->steps());

      return this->fieldMemory();
   }

   int ImExRK3::nonlinearMemory(const int step) const
   {
      assert(step < this->steps());

      int mem;

      // First step does not use t_(n-1) nonlinear term
      if(step == 0)
      {
         mem = 0;
      } else
      {
         mem = this->nonlinearMemory();
      }

      return mem;
   }

   void ImExRK3::init()
   {
   }

   bool ImExRK3::hasEmbedded() const
   {
      return false;
   }

}
}
