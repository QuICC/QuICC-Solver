/**
 * @file ISphericalInviscidCflWrapper.hpp
 * @brief CFL constraint in a spherical geometry
 */

#ifndef QUICC_DIAGNOSTICS_ISPHERICALINVISCIDCFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_ISPHERICALINVISCIDCFLWRAPPER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Diagnostics/ICflWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a spherical geometry
    */
   class ISphericalInviscidCflWrapper: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          */
         ISphericalInviscidCflWrapper(const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalInviscidCflWrapper() = default;

      protected:
         /**
          * @brief Get effective max harmonic degree L
          */
         virtual MHDFloat effectiveMaxL(const MHDFloat r) const = 0;

         /**
          * @brief Initialise the mesh spacings
          */
         void initMesh(const std::vector<Array>& mesh);

         /**
          * @brief Spacing between grid points
          */
         std::vector<Array> mMeshSpacings;

      private:
   };

   /// Typedef for a shared ISphericalInviscidCflWrapper
   typedef std::shared_ptr<ISphericalInviscidCflWrapper> SharedISphericalInviscidCflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_ISPHERICALINVISCIDCFLWRAPPER_HPP
