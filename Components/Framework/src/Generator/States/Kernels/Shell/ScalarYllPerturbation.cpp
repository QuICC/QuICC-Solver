/**
 * @file ScalarYllPerturbation.cpp
 * @brief Source of for scalar field perturbation with Y_l^l spherical harmonic generator kernel
 */

// System includes
//
#include <cmath>

// Project includes
//
#include "Types/Math.hpp"
#include "QuICC/Generator/States/Kernels/Shell/ScalarYllPerturbation.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Shell {

   ScalarYllPerturbation::ScalarYllPerturbation()
      : IPhysicalKernel()
   {
   }

   void ScalarYllPerturbation::init(const MHDFloat ri, const MHDFloat ro, const MHDFloat epsilon, const int l)
   {
      this->mRi = ri;
      this->mRo = ro;
      this->mEpsilon = epsilon;
      this->mL = l;
   }

   void ScalarYllPerturbation::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);

      MHDFloat ri = this->mRi;
      MHDFloat ro = this->mRo;
      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);

      const int& l = this->mL;
      const MHDFloat dl = static_cast<MHDFloat>(l);

      MHDFloat funcR;
      MHDFloat funcTh;

      // Normalize spherical harmonic Yll
      MHDFloat cYll = 1.0/(std::pow(2.0,l)*std::sqrt(4.0*Math::PI))*std::exp(0.5*std::lgamma((2*l+1) + 1) - std::lgamma(l + 1));
      const MHDFloat eps = this->mEpsilon*cYll;

      Array funcPh = (dl*phGrid).array().cos();

      MHDFloat r;
      MHDFloat x;
      MHDFloat theta;

      auto& tRes = *this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D);

      rNLComp.rData().setConstant(0);
      nR = tRes.dim<Dimensions::Data::DAT3D>();
      for(int iR = 0; iR < nR; ++iR)
      {
         r = rGrid(tRes.idx<Dimensions::Data::DAT3D>(iR));
         x = 2.0*r - ri - ro;
         funcR = 1.0 - 3.0*std::pow(x,2) + 3.0*std::pow(x,4) - std::pow(x,6);

         nTh = tRes.dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(tRes.idx<Dimensions::Data::DAT2D>(iTh, iR));
            funcTh = std::pow(std::sin(theta),l);

            rNLComp.addProfile(eps*funcR*funcTh*funcPh,iTh,iR);
         }
      }
   }

} // Shell
} // Kernel
} // Physical
} // QuICC
