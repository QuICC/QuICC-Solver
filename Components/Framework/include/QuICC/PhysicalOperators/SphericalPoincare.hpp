/*
 * @file SphericalPoincare.hpp
 * @brief Implementation of the spherical poincare term
 */

#ifndef QUICC_PHYSICAL_SPHERICALPOINCARE_HPP
#define QUICC_PHYSICAL_SPHERICALPOINCARE_HPP

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

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical poincare term
    */
   class SphericalPoincare
   {
      public:
         /**
          * @brief Set S to Poincare term
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Add Poincare term to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

         /**
          * @brief Substract Poincare term from S
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalPoincare() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalPoincare() = default;
   };
}
}

#endif // QUICC_PHYSICAL_SPHERICALPOINCARE_HPP
