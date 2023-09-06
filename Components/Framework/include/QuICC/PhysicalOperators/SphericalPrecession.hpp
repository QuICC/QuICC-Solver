/*
 * @file SphericalPrecession.hpp
 * @brief Implementation of the spherical Coriolis + precession term
 */

#ifndef QUICC_PHYSICAL_SPHERICALPRECESSION_HPP
#define QUICC_PHYSICAL_SPHERICALPRECESSION_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"

namespace QuICC {

namespace Physical {

   /**
    * @brief Implementation of the spherical Coriolis + precession term
    */
   class SphericalPrecession
   {
      public:
         /**
          * @brief Set S to Coriolis + precession term
          */
         static void set(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

         /**
          * @brief Add Coriolis + precesion term to S
          */
         static void add(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

         /**
          * @brief Substract (Coriolis + precession) term from S
          */
         static void sub(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c = 1.0);

      protected:

      private:
         /**
          * @brief Empty constructor
          */
         SphericalPrecession() = default;

         /**
          * @brief Empty destructor
          */
         ~SphericalPrecession() = default;
   };
}
}

#endif // QUICC_PHYSICAL_SPHERICALPRECESSION_HPP
