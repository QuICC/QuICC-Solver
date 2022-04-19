/**
 * @file ShellCflWrapper.hpp
 * @brief CFL constraint in a spherical shell
 */

#ifndef QUICC_DIAGNOSTICS_SHELLCFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_SHELLCFLWRAPPER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Diagnostics/ISphericalCflWrapper.hpp"
#include "QuICC/Diagnostics/IVectorWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a spherical Shell
    */
   class ShellCflWrapper: public ISphericalCflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          */
         ShellCflWrapper(const SharedIVectorWrapper spVelocity, const std::map<std::size_t,NonDimensional::SharedINumber>& params);

         /**
          * @brief Constructor
          *
          * @param spVelocity Velocity wrapper
          * @param spMagnetic Magnetic wrapper
          */
         ShellCflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic, const std::map<std::size_t,NonDimensional::SharedINumber>& params);

         /**
          * @brief Destructor
          */
         virtual ~ShellCflWrapper();

      protected:

      private:
         /**
          * @brief Get effective max harmonic degree L
          */
         virtual MHDFloat effectiveMaxL(const MHDFloat r) const;
   };

   /// Typedef for a shared ShellCflWrapper
   typedef std::shared_ptr<ShellCflWrapper> SharedShellCflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_SHELLCFLWRAPPER_HPP
