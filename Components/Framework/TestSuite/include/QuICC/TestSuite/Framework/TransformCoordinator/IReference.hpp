/**
 * @file TestHelper.hpp
 * @brief Helper functions to setup TransformCoordinator tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TCOORD_IREFERENCE_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TCOORD_IREFERENCE_HPP

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/Test.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   class IReference
   {
      public:
         // Constructor
         IReference() = default;

         // Destructor
         virtual ~IReference() = default;

         /// @brief Input for scalar field
         virtual MHDComplex inScalar(Test& test, int i, int j, int k) = 0;

         /// @brief Input for toroidal field
         virtual MHDComplex inTor(Test& test, int i, int j, int k) = 0;

         /// @brief Input for poloidal field
         virtual MHDComplex inPol(Test& test, int i, int j, int k) = 0;

         /// @brief Reference output for scalar field
         virtual MHDComplex refScalar(Test& test, int i, int j, int k) = 0;

         /// @brief Reference output for toroidal field
         virtual MHDComplex refTor(Test& test, int i, int j, int k) = 0;

         /// @brief Reference output for poloidal field
         virtual MHDComplex refPol(Test& test, int i, int j, int k) = 0;
   };

   /**
    * @brief Generate unit spectrum reference
    */
   MHDComplex unitReference(const Test& test, const int i, const int j, const int k);

   /**
    * @brief Generate unit spectrum reference for spherical harmonics
    */
   MHDComplex unitReferenceSH(const int n, const int l, const int m);

   /**
    * @brief Generate unit spectrum reference for double Fourier series
    */
   MHDComplex unitReferenceFF(const int n, const int k1, const int k2);
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TCOORD_IREFERENCE_HPP
