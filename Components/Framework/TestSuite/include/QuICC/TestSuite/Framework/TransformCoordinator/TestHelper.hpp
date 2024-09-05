/**
 * @file TestHelper.hpp
 * @brief Helper functions to setup TransformCoordinator tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TCOORD_TESTHELPER_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TCOORD_TESTHELPER_HPP

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/Communicators/Communicator.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingTools.hpp"
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingDescription.hpp"
#include "QuICC/Enums/VectorFormulation.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Enums/SplittingTools.hpp"
#include "QuICC/LoadSplitter/LoadSplitter.hpp"
#include "QuICC/Transform/Setup/Default.hpp"
#include "QuICC/TestSuite/Framework/TransformCoordinator/Test.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   // Typedef for storing result of error checks
   typedef std::tuple<bool,MHDFloat,MHDFloat> ErrorType;

   /**
    * @brief Create distributed resolution object
    *
    * @param dim        Spectral dimensions
    * @param algorithm  Splitting algorithm
    * @param grouper    Communication grouping algorithm
    * @param factors    Imposed CPU factorization
    * @param opt1D      Truncation ID for 1D
    */
   template <typename TScheme> std::pair<SharedResolution,Parallel::SplittingDescription> initResolution(ArrayI& dim, const std::string algorithm, const std::string grouper, const std::list<int>& factors, const std::vector<std::size_t>& opt1D);

   /**
    * @brief Process command line and generate dimension array
    */
   ArrayI processCmdLine(Test& test);

   /**
    * @brief Init variables
    */
   void initVariables(Test& test);

   /**
    * @brief Init physical computation kernels
    */
   void initKernels(Test& test);

   /**
    * @brief Initialize transform trees
    */
   void initTrees(Test& test);

   /**
    * @brief Initialize the transform coordinator
    */
   void initCoordinator(Test& test, const Parallel::SplittingDescription& descr);

   /**
    * @brief Set variables
    */
   void setVariables(Test& test);

   /**
    * @brief Set variables to bad value
    */
   void scrambleVariables(Test& test);

   /**
    * @brief Check variables
    */
   void checkVariables(Test& test);

   /**
    * @brief Backward transform from spectral to physical
    */
   void backward(Test& test);

   /**
    * @brief Nonlinear computation on physical grid and forward transform from physical to spectral
    */
   void nonlinearAndForward(Test& test);

   /**
    * @brief Compute ULP
    */
   ErrorType computeUlp(const MHDComplex data, const MHDComplex ref, MHDFloat refMod, const MHDFloat tol, const MHDFloat eps);

   /**
    * @brief Compute ULP
    */
   ErrorType computeUlp(const MHDFloat data, const MHDFloat ref, const MHDFloat refMod, const MHDFloat tol, const MHDFloat eps);

   template <typename TScheme> std::pair<SharedResolution,Parallel::SplittingDescription> initResolution(ArrayI& dim, const std::string algorithm, const std::string grouper, const std::list<int>& factors, const std::vector<std::size_t>& opt1D)
   {
      INFO( "MPI rank: " << QuICC::QuICCEnv().id() );
      INFO( "MPI size: " << QuICC::QuICCEnv().size() );
      QuICC::Parallel::LoadSplitter splitter(QuICC::QuICCEnv().id(), QuICC::QuICCEnv().size());

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
      auto best = splitter.bestSplitting(true);
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
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TCOORD_TESTHELPER_HPP
