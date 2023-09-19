/**
 * @file SphericalHarmonic.hpp
 * @brief Class to hold tools for spherical harmonic computations in physical space
 */

#ifndef QUICC_EQUATIONS_SPHERICALHARMONIC_HPP
#define QUICC_EQUATIONS_SPHERICALHARMONIC_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Tools {

   /**
    * @brief Class to hold the list of possible exact states
    */
   struct SphericalHarmonic
   {
      /**
       * @brief Compute spherical harmonic physical values
       */
      static void Ylm(Array& rField, const int l, const int m, const MHDComplex c, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Toroidal Y_l^m
       */
      static void Torlm(Array& rField, const int l, const int m, FieldComponents::Physical::Id compId, const MHDComplex ct, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Poloidal Y_l^m
       */
      static void Pollm(Array& rField, const int l, const int m, FieldComponents::Physical::Id compId, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Toroidal Y_0^0
       */
      static void Tor00(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex ct, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Toroidal Y_1^0
       */
      static void Tor10(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex ct, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Toroidal Y_1^1
       */
      static void Tor11(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex ct, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Toroidal Y_2^0
       */
      static void Tor20(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex ct, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Toroidal Y_2^1
       */
      static void Tor21(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex ct, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Toroidal Y_2^2
       */
      static void Tor22(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex ct, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Toroidal Y_5^4
       */
      static void Tor54(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex ct, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Poloidal Y_0^0
       */
      static void Pol00(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Poloidal Y_1^0
       */
      static void Pol10(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Poloidal Y_1^1
       */
      static void Pol11(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Poloidal Y_2^0
       */
      static void Pol20(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Poloidal Y_2^1
       */
      static void Pol21(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Poloidal Y_2^2
       */
      static void Pol22(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi);

      /**
       * @brief Generate Poloidal Y_4^3
       */
      static void Pol43(Array& rField, FieldComponents::Physical::Id compId, const MHDComplex cq, const MHDComplex cs, const MHDFloat theta, const Array& phi);

   };

}
}
}
}

#endif // QUICC_EQUATIONS_SPHERICALHARMONIC_HPP
