/**
 * @file TubularSplitting.hpp
 * @brief Implementation of a double load splitting algorithm, aka "Tubular" splitting
 */

#ifndef QUICC_PARALLEL_TUBULARSPLITTING_HPP
#define QUICC_PARALLEL_TUBULARSPLITTING_HPP

// System includes
//
#include <utility>
#include <vector>

// Project includes
//
#include "QuICC/Enums/Splitting.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "QuICC/Resolutions/TransformResolution.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of a double load splitting algorithm, aka "Tubular" splitting
    */
   class TubularSplitting: public SplittingAlgorithm
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id ID of the CPU/Core
          * @param nCpu Number of cores used
          * @param dim  Dimensions
          * @param algorithm  Splitting algorithm
          * @param factors Imposed CPU factorizations
          */
         TubularSplitting(const int id, const int nCpu, const ArrayI& dim, Splitting::Algorithms::Id algorithm, const std::list<int>& factors);

         /**
          * @brief Destructor
          */
         ~TubularSplitting() = default;

         /**
          * @brief Check if factorisation is applicable to scheme
          */
         bool applicable() const final;

      protected:
         /**
          * @brief Split ith dimension transform
          *
          * @param transId Split the ith dimension
          * @param cpuId   ID of the CPU
          * @param status  Status output
          */
         SharedTransformResolution splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status) final;

         /**
          * @brief Select the transform grouper
          */
         void selectGrouper(const Splitting::Groupers::Id selected);

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         Array computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp) final;

      private:
   };

}
}

#endif // QUICC_PARALLEL_TUBULARSPLITTING_HPP
