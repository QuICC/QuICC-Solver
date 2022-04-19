/**
 * @file SphericalTorPolWrapper.hpp
 * @brief Spherical Toroidal/Poloidal decomposition wrapper into velocity field
 */

#ifndef QUICC_DIAGNOSTICS_SPHERICALTORPOLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_SPHERICALTORPOLWRAPPER_HPP

// Configuration includes
//

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Diagnostics/IVectorWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Spherical Toroidal/Poloidal decomposition wrapper into velocity field
    */
   class SphericalTorPolWrapper: public IVectorWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         SphericalTorPolWrapper(const Framework::Selector::VariantSharedVectorVariable spTorPol);

         /**
          * @brief Destructor
          */
         virtual ~SphericalTorPolWrapper();

         /**
          * @brief Get first velocity field component
          */
         virtual const Framework::Selector::PhysicalScalarField& one() const;

         /**
          * @brief Get second velocity field component
          */
         virtual const Framework::Selector::PhysicalScalarField& two() const;

         /**
          * @brief Get third velocity field component
          */
         virtual const Framework::Selector::PhysicalScalarField& three() const;

         /**
          * @brief Get Resolution
          */
         virtual const Resolution& res() const;

      protected:

      private:
         /**
          * @brief Shared velocity vector variable
          */
         Framework::Selector::VariantSharedVectorVariable mspTorPol;
   };

   /// Typedef for a shared SphericalTorPolWrapper
   typedef std::shared_ptr<SphericalTorPolWrapper> SharedSphericalTorPolWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_SPHERICALTORPOLWRAPPER_HPP
