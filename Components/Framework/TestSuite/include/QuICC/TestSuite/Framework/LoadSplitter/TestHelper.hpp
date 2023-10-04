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
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/Enums/VectorFormulation.hpp"
#include "QuICC/Transform/Setup/Default.hpp"
#include "QuICC/Transform/Setup/GaussianQuadrature.hpp"
#include "QuICC/Transform/Setup/Uniform.hpp"
#include "QuICC/Transform/Setup/Triangular.hpp"
#include "QuICC/Transform/Setup/Trapezoidal.hpp"
#include "QuICC/Enums/SplittingTools.hpp"
#include "QuICC/LoadSplitter/LoadSplitter.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
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
    * @param sRes       Simulation resolution
    * @param tRes       Transform resolution
    */
   void writeData(const std::string fname, const SimRes& sRes, const TransformResolution& tRes);

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
    * @param sRes    Simulation resolution
    * @param tRes    Transform resolution
    */
   void checkReference(const std::string fname, const SimRes& sRes, const TransformResolution& tRes);

   /**
    * @brief Check transform resolution against full serial reference data
    *
    * @param fname   Filename of reference
    * @param tRes    Transform resolution
    */
   void checkSerialReference(const std::string fname, const std::multimap<int,int>& modes);

   /**
    * @brief Create distributed resolution object
    *
    * @param dim        Spectral dimensions
    * @param algorithm  Splitting algorithm
    * @param grouper    Communication grouping algorithm
    * @param factors    Imposed CPU factorization
    * @param opt1D      Truncation ID for 1D
    */
   template <typename TScheme> std::pair<SharedResolution,Parallel::SplittingDescription> initResolution(const int id, const int np, ArrayI& dim, const std::string algorithm, const std::string grouper, const std::list<int>& factors, const std::vector<std::size_t>& opt1D);


   template <typename TScheme> std::pair<SharedResolution,Parallel::SplittingDescription> initResolution(const int id, const int np, ArrayI& dim, const std::string algorithm, const std::string grouper, const std::list<int>& factors, const std::vector<std::size_t>& opt1D)
   {
      QuICC::Parallel::LoadSplitter splitter(id, np);

      // Create spatial scheme
      auto spScheme = std::make_shared<TScheme>(QuICC::VectorFormulation::TORPOL, QuICC::GridPurpose::SIMULATION);
      std::map<std::size_t,std::vector<std::size_t>> transformSetup;

      // Set 1D implementation
      transformSetup.emplace(0, opt1D);
      // Set other dimensions
      for(int i = 1; i < dim.size(); i++)
      {
         std::vector<std::size_t> opt = {QuICC::Transform::Setup::Default::id()};
         transformSetup.emplace(i, opt);
      }
      spScheme->setImplementation(transformSetup);

      // Create scheme's builder
      auto spBuilder = spScheme->createBuilder(dim, false);

      // Select splitting algorithm
      auto algoId = QuICC::Splitting::getAlgorithmId(algorithm);
      auto grouperId = QuICC::Splitting::getGrouperId(grouper);
      INFO( "Comm algorithm: " << algorithm );
      INFO( "Comm grouper: " << grouper );
      std::set<QuICC::Splitting::Algorithms::Id> algos = {algoId};
      splitter.init(spBuilder, algos, grouperId, factors);

      // Generate resolution
      auto best = splitter.bestSplitting(false);
      auto spRes = best.first;
      spBuilder->tuneResolution(spRes, best.second);

      // Set boxscale
      Array box(3);
      box << 1.0,1.0,1.0;
      spRes->setBoxScale(box);

      // Pass spatial scheme to resolution
      spRes->setSpatialScheme(spScheme);

      return best;
   }

} // LoadSplitter
} // Framework
} // TestSuite
} // QuIIC

#endif //QUICC_TESTSUITE_FRAMEWORK_LOADSPLITTER_TESTHELPER_HPP
