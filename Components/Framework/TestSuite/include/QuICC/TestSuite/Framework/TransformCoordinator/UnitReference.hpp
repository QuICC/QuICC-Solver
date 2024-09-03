/**
 * @file TestHelper.hpp
 * @brief Helper functions to setup TransformCoordinator tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TCOORD_UNITREFERENCE_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TCOORD_UNITREFERENCE_HPP

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/IReference.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/Test.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   class UnitReference: public IReference
   {
      public:
         // Constructor
         UnitReference(const Test::FieldId id);

         // Destructor
         ~UnitReference() = default;

         /// @brief Input for scalar field
         MHDComplex inScalar(Test& test, int i, int j, int k) final;

         /// @brief Input for toroidal field
         MHDComplex inTor(Test& test, int i, int j, int k) final;

         /// @brief Input for poloidal field
         MHDComplex inPol(Test& test, int i, int j, int k) final;

         /// @brief Reference output for scalar field
         MHDComplex refScalar(Test& test, int i, int j, int k) final;

         /// @brief Reference output for toroidal field
         MHDComplex refTor(Test& test, int i, int j, int k) final;

         /// @brief Reference output for poloidal field
         MHDComplex refPol(Test& test, int i, int j, int k) final;

      private:
         Test::FieldId mFieldId;
   };
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TCOORD_UNITREFERENCE_HPP
