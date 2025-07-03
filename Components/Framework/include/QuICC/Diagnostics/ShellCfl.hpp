/**
 * @file ShellCfl.hpp
 * @brief CFL constraint in a spherical shell
 */

#ifndef QUICC_DIAGNOSTICS_SHELLCFL_HPP
#define QUICC_DIAGNOSTICS_SHELLCFL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Diagnostics/ISphericalCflWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a spherical Shell
    */
   template <typename TCfl>
   class ShellCfl: public TCfl
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          */
         ShellCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~ShellCfl() = default;

      protected:

      private:
         /**
          * @brief Get effective max harmonic degree L
          */
         virtual MHDFloat effectiveMaxL(const MHDFloat r) const;
   };

   template <typename TCfl>
   ShellCfl<TCfl>::ShellCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant)
      : TCfl(params, courant)
   {
   }

   template <typename TCfl>
   MHDFloat ShellCfl<TCfl>::effectiveMaxL(const MHDFloat r) const
   {
      const auto& field = this->mFields.begin()->second;
      return static_cast<MHDFloat>(field->res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL) - 1);
   }
}
}

#endif // QUICC_DIAGNOSTICS_SHELLCFL_HPP
