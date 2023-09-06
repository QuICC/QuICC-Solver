/**
 * @file TestHelper.hpp
 * @brief Helper functions for test
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTHELPER_HPP
#define QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTHELPER_HPP

// System includes
//
#include <vector>
#include <string>

// Project includes
//
#include "QuICC/Resolutions/TransformResolution.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace LoadSplitter {

   /**
    * Data structure to pass simulation wide truncation
    */
   struct SimRes
   {
      /// number of spectral modes
      int nSpec;
      /// number of physical grid points
      int nPhys;
   };

   /**
    * @brief Write resolution data to file
    *
    * @param fname      Filename
    * @param tRes       Resolution
    * @param isDetailed Write detailed output? (summary vs mode list)
    */
   void writeData(const std::string fname, const SimRes& sRes, const TransformResolution& tRes, const bool isDetailed);

   /**
    * @brief Read reference data
    *
    * @param path File path
    * @param data Storage for read data
    */
   void readData(const std::string path, std::vector<int>& data);

   /**
    * @brief Collect (i,j) modes
    *
    * @param tRes    Transform resolution
    * @param modes   Storage for modes
    */
   void collectModes(const TransformResolution& tRes, std::multimap<int,int>& modes);

   /**
    * @brief Check transform resolution against distributed reference data
    *
    * @param fname   Filename of reference
    * @param tRes    Transform resolution
    */
   void checkReference(const std::string fname, const TransformResolution& tRes);

   /**
    * @brief Check transform resolution against full serial reference data
    *
    * @param fname   Filename of reference
    * @param tRes    Transform resolution
    */
   void checkSerialReference(const std::string fname, const std::multimap<int,int>& modes);

} // LoadSplitter
} // Framework
} // TestSuite
} // QuIIC

#endif //QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTHELPER_HPP
