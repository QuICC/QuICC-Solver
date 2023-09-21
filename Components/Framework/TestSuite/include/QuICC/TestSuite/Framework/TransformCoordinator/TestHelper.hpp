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
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingDescription.hpp"
#include "QuICC/Enums/VectorFormulation.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Enums/SplittingTools.hpp"
#include "QuICC/LoadSplitter/LoadSplitter.hpp"
#include "QuICC/Communicators/Communicator.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TransformCoordinators/TransformCoordinator.hpp"
#include "QuICC/TransformGroupers/IForwardGrouper.hpp"
#include "QuICC/TransformGroupers/IBackwardGrouper.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/VectorFields/VectorField.hpp"
#include "QuICC/Transform/Setup/Default.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   // Typedef for storing result of error checks
   typedef std::tuple<bool,MHDFloat,MHDFloat> ErrorType;

   /**
    * @brief Small struct to collect the different objects needed
    */
   struct Test
   {
      enum class FieldId
      {
         SCALAR= 0,
         TOR,
         POL,
         TORPOL,
         SCALAR_AND_TORPOL,
      };

      enum class KernelId
      {
         PASSTHROUGH = 0,
      };

      enum class PathId
      {
         BFLOOP = 0,
      };

      enum class SpectrumId
      {
         UNIT = 0,
      };

      /**
       * @brief construtor
       */
      Test();

      /**
       * @brief Translate ID to test configuration
       */
      void configure(const int id);

      /**
       * @brief Tolerance for error checks
       */
      MHDFloat tolerance() const;

      /**
       * @brief Shared resolution
       */
      SharedResolution spRes;

      /**
       * @brief TransformCoordinator
       */
      TransformCoordinator<Parallel::Communicator> coord;

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

      /**
       * @brief Scalar variables
       */
      std::map<std::size_t, ::QuICC::Framework::Selector::VariantSharedScalarVariable> scalars;

      /**
       * @brief Vector variables
       */
      std::map<std::size_t, ::QuICC::Framework::Selector::VariantSharedVectorVariable> vectors;

      /**
       * @brief Physical space computational kernels
       */
      std::map<std::size_t, Physical::Kernel::SharedIPhysicalKernel> kernels;

      /**
       * @brief Epsilon for error checks
       */
      MHDFloat epsilon;

      /**
       * @brief Max acceptable ULP for error checks
       */
      MHDFloat maxUlp;

      /**
       * @brief Test field ID
       */
      FieldId fieldId;

      /**
       * @brief Test kernel ID
       */
      KernelId kernelId;

      /**
       * @brief Test transform path ID
       */
      PathId pathId;

      /**
       * @brief Test input spectrum ID
       */
      SpectrumId spectrumId;
   };

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
   ArrayI processCmdLine();

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

      Array box(3);
      box << 1.0,1.0,1.0;
      spRes->setBoxScale(box);

      spRes->setSpatialScheme(spScheme);

      return best;
   }
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TCOORD_TESTHELPER_HPP
