/**
 * @file SphereCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a full sphere
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/SphereCflWrapper.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

namespace QuICC {

namespace Diagnostics {

   SphereCflWrapper::SphereCflWrapper(const SharedIVectorWrapper spVelocity, const std::map<std::size_t,NonDimensional::SharedINumber>& params)
      : ISphericalCflWrapper(spVelocity, params)
   {
   }

   SphereCflWrapper::SphereCflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic, const std::map<std::size_t,NonDimensional::SharedINumber>& params)
      : ISphericalCflWrapper(spVelocity, spMagnetic, params)
   {
   }

   MHDFloat SphereCflWrapper::effectiveMaxL(const MHDFloat r) const
   {
      MHDFloat l = 0;
      for(; l < static_cast<MHDFloat>(this->mspVelocity->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)); l++)
      {
         if(r < this->jacobiRoot(l))
         {
            break;
         }
      }

      return l;
   }

   MHDFloat SphereCflWrapper::jacobiRoot(const MHDFloat l) const
   {
      MHDFloat rN = std::ceil(3.*this->mspVelocity->res().counter().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, l)/2.);

      MHDFloat a = -0.5;
      MHDFloat b = l - 0.5;
      MHDFloat rho = (rN + l/2.0);

      // Compute 4th order approximation to Bessel zero
      MHDFloat jb = b + 1.8557571*std::cbrt(b);
      jb += 1.033150*(1./std::cbrt(b));
      jb -= 0.00397*std::pow(b,-1.);
      jb -= 0.0908*std::pow(std::cbrt(b),-5.);
      jb += 0.043*std::pow(std::cbrt(b),-7.);

      // Compute 2nd order approximation to first Jacobi root
      MHDFloat phi = jb/rho;
      MHDFloat theta = phi;
      theta += ((b*b - 0.25)*(1. - phi*(1./std::tan(phi)))/(2.*phi) - 0.25*(b*b - a*a)*std::tan(0.5*phi))/(rho*rho);

      return std::sin(theta/2.0);
   }

}
}
