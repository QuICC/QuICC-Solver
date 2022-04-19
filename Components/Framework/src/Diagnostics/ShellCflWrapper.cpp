/**
 * @file ShellCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a spherical shell
 */

// Debug includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Diagnostics/ShellCflWrapper.hpp"

// Project includes
//

namespace QuICC {

namespace Diagnostics {

   ShellCflWrapper::ShellCflWrapper(const SharedIVectorWrapper spVelocity, const std::map<std::size_t,NonDimensional::SharedINumber>& params)
      : ISphericalCflWrapper(spVelocity, params)
   {
   }

   ShellCflWrapper::ShellCflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic, const std::map<std::size_t,NonDimensional::SharedINumber>& params)
      : ISphericalCflWrapper(spVelocity, spMagnetic, params)
   {
   }

   ShellCflWrapper::~ShellCflWrapper()
   {
   }

   MHDFloat ShellCflWrapper::effectiveMaxL(const MHDFloat r) const
   {
      return static_cast<MHDFloat>(this->mspVelocity->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL) - 1);
   }

}
}
