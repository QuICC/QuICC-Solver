/**
 * @file TestHelper.hpp
 * @brief Helper functions to setup TransformCoordinator tests
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_TCOORD_TEST_HPP
#define QUICC_TESTSUITE_FRAMEWORK_TCOORD_TEST_HPP

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

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace TCoord {

   // Typedef for storing result of error checks
   typedef std::tuple<bool,MHDFloat,MHDFloat> ErrorType;

   class IReference;

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
         INERTIA,
         CURL
      };

      enum class PathId
      {
         BFLOOP = 0,
         NLLOOP,
      };

      enum class SpectrumId
      {
         UNIT = 0,
         SINGLE_MODE,
         FROM_FILE
      };

      /**
       * @brief construtor
       */
      Test();

      /**
       * @brief Translate ID to test configuration
       *
       * @param id      Test ID
       * @param subdir  Reference data subdirectory
       */
      void configure(const int id, const std::string subdir);

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

      std::shared_ptr<IReference> spRef;

      /**
       * @brief Filename base
       */
      std::string fbase;
   };
}
}
}
}

#endif //QUICC_TESTSUITE_FRAMEWORK_TCOORD_TEST_HPP
