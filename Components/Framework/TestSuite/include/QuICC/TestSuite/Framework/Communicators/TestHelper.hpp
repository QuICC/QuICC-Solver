/**
 * @file TestHelper.hpp
 * @brief Different helper functions for Communicator test
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_COMMUNICATORS_TESTHELPER_HPP
#define QUICC_TESTSUITE_FRAMEWORK_COMMUNICATORS_TESTHELPER_HPP

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Enums/SplittingTools.hpp"
#include "QuICC/LoadSplitter/LoadSplitter.hpp"
#include "QuICC/Communicators/Communicator.hpp"
#include "QuICC/TransformGroupers/IForwardGrouper.hpp"
#include "QuICC/TransformGroupers/IBackwardGrouper.hpp"
#include "QuICC/Transform/Setup/Default.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Communicators {

   /**
    * @brief Small struct to collect the different objects needed to perform all transpose stages
    */
   struct Test
   {
      /**
       * @brief Shared resolution
       */
      SharedResolution spRes;

      /**
       * @brief Communicator
       */
      Parallel::Communicator comm;

      /**
       * @brief Forward grouper
       */
      Transform::SharedIForwardGrouper spFwdGrouper;

      /**
       * @brief Backward grouper
       */
      Transform::SharedIBackwardGrouper spBwdGrouper;

      /**
       * @brief Backward transform tree
       */
      std::vector<Transform::TransformTree> bwdTree;

      /**
       * @brief Forward transform tree
       */
      std::vector<Transform::TransformTree> fwdTree;
   };

   /**
    * @brief Constant for unused value for output data
    */
   static const MHDFloat badValueOut = std::numeric_limits<QuICC::MHDFloat>::max();

   /**
    * @brief Constant for unused value for input data
    */
   static const MHDFloat badValueIn = std::numeric_limits<QuICC::MHDFloat>::max()/2.0;

   /**
    * @brief Create distributed resolution object
    *
    * @param dim        Spectral dimensions
    * @param algorithm  Splitting algorithm
    * @param grouper    Communication grouping algorithm
    * @param factors    Imposed CPU factorization
    */
   template <typename TScheme> std::pair<SharedResolution,Parallel::SplittingDescription> initResolution(ArrayI& dim, const std::string algorithm, const std::string grouper, const std::list<int>& factors);

   /**
    * @brief Process command line options
    */
   ArrayI processCmdLine();

   /**
    * @brief Initialize the communicator
    */
   void initCommunicator(Test& test, const Parallel::SplittingDescription& descr);

   /**
    * @brief Setup commmunication
    */
   void setupSpectralCommunication(Test& test);

   /**
    * @brief Perform transform of data between 1D and 2D transform
    */
   void transposeSpectral_1D(Test& test);

   /**
    * @brief Setup commmunication
    */
   void setup1D2D3DCommunication(Test& test);

   /**
    * @brief Perform transform of data between 1D and 2D transform
    */
   void transpose1D_2D(Test& test);

   /**
    * @brief Perform transform of data between 2D and 3D transform
    */
   void transpose2D_3D(Test& test);

   template <typename TScheme> std::pair<SharedResolution,Parallel::SplittingDescription> initResolution(ArrayI& dim, const std::string algorithm, const std::string grouper, const std::list<int>& factors)
   {
      INFO( "MPI rank: " << QuICC::QuICCEnv().id() );
      INFO( "MPI size: " << QuICC::QuICCEnv().size() );
      QuICC::Parallel::LoadSplitter splitter(QuICC::QuICCEnv().id(), QuICC::QuICCEnv().size());

      // Create spatial scheme
      auto spScheme = std::make_shared<TScheme>(QuICC::VectorFormulation::TORPOL, QuICC::GridPurpose::SIMULATION);
      std::map<std::size_t,std::vector<std::size_t>> transformSetup;
      for(int i = 0; i < dim.size(); i++)
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
      auto best = splitter.bestSplitting(true);
      auto spRes = best.first;
      spBuilder->tuneResolution(spRes, best.second);
      spRes->setSpatialScheme(spScheme);

      return best;
   }
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_COMMUNICATORS_TESTHELPER_HPP
