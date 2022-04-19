/**
 * @file CartesianTorPolWrapper.hpp
 * @brief Cartesian Toroidal/Poloidal decomposition wrapper into velocity field
 */

#ifndef QUICC_DIAGNOSTICS_CARTESIANTORPOLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_CARTESIANTORPOLWRAPPER_HPP

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
    * @brief Cartesian Toroidal/Poloidal decomposition wrapper into velocity field
    */
   class CartesianTorPolWrapper: public IVectorWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         CartesianTorPolWrapper(const Framework::Selector::VariantSharedVectorVariable spTorPol);

         /**
          * @brief Destructor
          */
         virtual ~CartesianTorPolWrapper();

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

   /// Typedef for a shared CartesianTorPolWrapper
   typedef std::shared_ptr<CartesianTorPolWrapper> SharedCartesianTorPolWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_CARTESIANTORPOLWRAPPER_HPP
