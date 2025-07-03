/**
 * @file SphereCfl.hpp
 * @brief CFL constraint in a full sphere
 */

#ifndef QUICC_DIAGNOSTICS_SPHERECFL_HPP
#define QUICC_DIAGNOSTICS_SPHERECFL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a full sphere
    */
   template <typename TCfl>
   class SphereCfl: public TCfl
   {
      public:
         /**
          * @brief Constructor
          */
         SphereCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~SphereCfl() = default;

      protected:

      private:
         /**
          * @brief Get effective max harmonic degree L
          */
         virtual MHDFloat effectiveMaxL(const MHDFloat r) const;

         /**
          * @brief Compute approximation to first Jacobi root for CFL
          */
         MHDFloat jacobiRoot(const MHDFloat l) const;
   };

   template <typename TCfl>
   SphereCfl<TCfl>::SphereCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant)
      : TCfl(params, courant)
   {
   }

   template <typename TCfl>
   MHDFloat SphereCfl<TCfl>::effectiveMaxL(const MHDFloat r) const
   {
      const auto& field = this->mFields.begin()->second;

      MHDFloat l = 0;
      for(; l < static_cast<MHDFloat>(field->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)); l++)
      {
         if(r < this->jacobiRoot(l))
         {
            break;
         }
      }

      return l;
   }

   template <typename TCfl>
   MHDFloat SphereCfl<TCfl>::jacobiRoot(const MHDFloat l) const
   {
      const auto& field = this->mFields.begin()->second;
      MHDFloat rN = std::ceil(3.*field->res().counter().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, l)/2.);

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

#endif // QUICC_DIAGNOSTICS_SPHERECFL_HPP
