/**
 * @file FileReference.cpp
 * @brief Read reference input and ouput from file
 */

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <fstream>
#include <limits>

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/FileReference.hpp"
#include "TestSuite/Io.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   FileReference::FileReference(const std::string path)
      :mPath(path)
   {
   }

   void FileReference::readFile(MatrixZ& data, const std::string path, const Test& test)
   {
      //auto nN = test.spRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      auto nL = test.spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
      auto nM = test.spRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
      int nModes = 0;
      for(int m = 0; m < nM; m++)
      {
         nModes += nL-m;
      }

      data.resize(nM, nModes);
      readData(data, path);
   }

   MHDComplex FileReference::getValue(MatrixZ& data, const std::string comp, const Test& test, const int i, const int j, const int k)
   {
      auto l = k;
      auto m = j;
      if(data.size() == 0)
      {
         std::string path = this->mPath + comp + ".dat";
         this->readFile(data, path, test);
      }

      //auto nL = test.spRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
      //auto nM = test.spRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
      int col = 0;
      for(int il = 0; il < l; il++)
      {
         col += il+1;
      }
      col += m;

      return data(i, col);
   }

   MHDComplex FileReference::inScalar(Test& test, int i, int j, int k)
   {
      return this->getValue(this->mInScalar, "_scalar_in", test, i, j, k);
   }

   MHDComplex FileReference::inTor(Test& test, int i, int j, int k)
   {
      return this->getValue(this->mInTor, "_tor_in", test, i, j, k);
   }

   MHDComplex FileReference::inPol(Test& test, int i, int j, int k)
   {
      return this->getValue(this->mInPol, "_pol_in", test, i, j, k);
   }

   MHDComplex FileReference::refScalar(Test& test, int i, int j, int k)
   {
      return this->getValue(this->mRefScalar, "_scalar_ref", test, i, j, k);
   }

   MHDComplex FileReference::refTor(Test& test, int i, int j, int k)
   {
      return this->getValue(this->mRefTor, "_tor_ref", test, i, j, k);
   }

   MHDComplex FileReference::refPol(Test& test, int i, int j, int k)
   {
      return this->getValue(this->mRefPol, "_pol_ref", test, i, j, k);
   }
}
}
}
}
