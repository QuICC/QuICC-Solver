/**
 * @file FileReference.hpp
 * @brief Helper functions to setup TransformCoordinator tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TCOORD_FILEREFERENCE_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TCOORD_FILEREFERENCE_HPP

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/Test.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/IReference.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   class FileReference: public IReference
   {
      public:
         // Constructor
         FileReference(const std::string path);

         // Destructor
         ~FileReference() = default;

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
         /// Read data file
         void readFile(MatrixZ& data, const std::string path, const Test& test);
         MHDComplex getValue(MatrixZ& data, const std::string comp, const Test& test, const int i, const int j, const int k);

         std::string mPath;
         MatrixZ mInScalar;
         MatrixZ mInTor;
         MatrixZ mInPol;
         MatrixZ mRefScalar;
         MatrixZ mRefTor;
         MatrixZ mRefPol;
   };
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TCOORD_FILEREFERENCE_HPP
