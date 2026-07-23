/**
 * @file SphericalBuoyancyAnelastic.hpp
 * @brief Implementation of the spherical buoyancy term
 */

#ifndef QUICC_PHYSICAL_SPHERICALBUOYANCYANELASTIC_HPP
#define QUICC_PHYSICAL_SPHERICALBUOYANCYANELASTIC_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "DenseSM/IGenericProfile.hpp"


namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical buoyancy term
    */
   class SphericalBuoyancyAnelastic
   {
      public:

         /**
          * @brief Add to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, 
                         FieldComponents::Physical::Id compId, 
                         const Resolution& res, 
                         const Array& r, 
                         const Framework::Selector::PhysicalScalarField &q, 
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pAlpha,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pCp,
                         std::shared_ptr<QuICC::DenseSM::IGenericProfile> pG,
                         const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalBuoyancyAnelastic() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalBuoyancyAnelastic() = default;
   };
}
}

#endif // QUICC_PHYSICAL_SPHERICALBUOYANCYANELASTIC_HPP
