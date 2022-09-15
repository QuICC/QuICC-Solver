/**
 * @file CommunicatorStages.hpp
 * @brief Different stages for Communicator test
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_COMMUNICATORS_COMMUNICATORSTAGES_HPP
#define QUICC_TESTSUITE_FRAMEWORK_COMMUNICATORS_COMMUNICATORSTAGES_HPP

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

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Communicators {

   /**
    * @brief Small struct to collect the different objects needed to perform all transpose stages
    */
   struct CommCoordinator
   {
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
    */
   template <typename TScheme> std::pair<SharedResolution,Parallel::SplittingDescription> initResolution(ArrayI& dim, const std::string algorithm, const std::string grouper);

   /**
    * @brief Initialize the communicator
    */
   void initCommunicator(CommCoordinator& coord, SharedResolution spRes, const Parallel::SplittingDescription& descr);

   /**
    * @brief Perform transform of data between 1D and 2D transform
    */
   void transposeSpectral_1D(const QuICC::Resolution& res, QuICC::Parallel::Communicator& comm);

   /**
    * @brief Perform transform of data between 1D and 2D transform
    */
   void transpose1D_2D(const QuICC::Resolution& res, QuICC::Parallel::Communicator& comm);

   /**
    * @brief Perform transform of data between 2D and 3D transform
    */
   void transpose2D_3D(const QuICC::Resolution& res, QuICC::Parallel::Communicator& comm);

   template <typename TScheme> std::pair<SharedResolution,Parallel::SplittingDescription> initResolution(ArrayI& dim, const std::string algorithm, const std::string grouper)
   {
      INFO( "MPI rank: " << QuICC::QuICCEnv().id() );
      INFO( "MPI size: " << QuICC::QuICCEnv().size() );
      QuICC::Parallel::LoadSplitter splitter(QuICC::QuICCEnv().id(), QuICC::QuICCEnv().size());
      auto spScheme = std::make_shared<TScheme>(QuICC::VectorFormulation::TORPOL, QuICC::GridPurpose::SIMULATION);
      auto spBuilder = spScheme->createBuilder(dim, false);

      // Select splitting algorithm
      auto algoId = QuICC::Splitting::getAlgorithmId(algorithm);
      auto grouperId = QuICC::Splitting::getGrouperId(grouper);
      INFO( "Comm algorithm: " << algorithm );
      INFO( "Comm grouper: " << grouper );
      std::set<QuICC::Splitting::Algorithms::Id> algos = {algoId};
      splitter.init(spBuilder, algos, grouperId);

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

#endif //QUICC_TESTSUITE_FRAMEWORK_COMMUNICATORS_COMMUNICATORSTAGES_HPP
