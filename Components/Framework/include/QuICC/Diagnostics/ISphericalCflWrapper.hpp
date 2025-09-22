/**
 * @file ISphericalCflWrapper.hpp
 * @brief CFL constraint in a spherical geometry
 */

#ifndef QUICC_DIAGNOSTICS_ISPHERICALCFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_ISPHERICALCFLWRAPPER_HPP

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
   class ISphericalCflWrapper: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          */
         ISphericalCflWrapper(const MHDFloat courant);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalCflWrapper() = default;

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

   /// Typedef for a shared ISphericalCflWrapper
   typedef std::shared_ptr<ISphericalCflWrapper> SharedISphericalCflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_ISPHERICALCFLWRAPPER_HPP
