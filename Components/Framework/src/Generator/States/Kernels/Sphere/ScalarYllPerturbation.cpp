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
#include "QuICC/Generator/States/Kernels/Sphere/ScalarYllPerturbation.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   ScalarYllPerturbation::ScalarYllPerturbation()
      : IPhysicalKernel()
   {
   }

   void ScalarYllPerturbation::init(const MHDFloat amplitude_bg, const MHDFloat epsilon, const int l)
   {
      this->mBg = amplitude_bg;
      this->mEpsilon = epsilon;
      this->mL = l;
   }

   void ScalarYllPerturbation::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);

      const int& l = this->mL;
      const MHDFloat dl = static_cast<MHDFloat>(l);

      // Normalize spherical harmonic Yll
      MHDFloat cYll = 1.0/(std::pow(2.0,l)*std::sqrt(4.0*Math::PI))*std::exp(0.5*std::lgamma((2*l+1) + 1) - std::lgamma(l + 1));

      MHDFloat funcR0;
      MHDFloat funcR3;
      MHDFloat funcTh0;
      MHDFloat funcTh3;
      Array funcPh0 = Array::Ones(nPh);
      Array funcPh3 = (dl*phGrid).array().cos() + (dl*phGrid).array().sin();

      // Background state is not solved for (no source term but background is imposed explicitly)
      MHDFloat eps = cYll*this->mEpsilon;

      MHDFloat r;
      MHDFloat theta;

      auto& tRes = *this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D);

      rNLComp.rData().setConstant(0);
      nR = tRes.dim<Dimensions::Data::DAT3D>();
      for(int iR = 0; iR < nR; ++iR)
      {
         r = rGrid(tRes.idx<Dimensions::Data::DAT3D>(iR));
         funcR0 = (1.0 - std::pow(r,2));
         funcR3 = std::pow(r,l)*(1.0 - std::pow(r,2));

         nTh = tRes.dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(tRes.idx<Dimensions::Data::DAT2D>(iTh, iR));
            funcTh0 = 1.0;
            funcTh3 = std::pow(std::sin(theta),l);

            // background state
            rNLComp.addProfile(this->mBg*funcR0*funcTh0*funcPh0,iTh,iR);

            // perturbation
            rNLComp.addProfile(eps*funcR3*funcTh3*funcPh3,iTh,iR);
         }
      }
   }

} // Sphere
} // Kernel
} // Physical
} // QuICC
