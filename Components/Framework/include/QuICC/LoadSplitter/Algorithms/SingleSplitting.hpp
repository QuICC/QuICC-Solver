/** 
 * @file SingleSplitting.hpp
 * @brief Implementation of a single load splitting algorithm
 */

#ifndef QUICC_PARALLEL_SINGLESPLITTING_HPP
#define QUICC_PARALLEL_SINGLESPLITTING_HPP

// System includes
//
#include <utility>
#include <vector>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "QuICC/Resolutions/TransformResolution.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of a single load splitting algorithm
    */
   class SingleSplitting: public SplittingAlgorithm
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id      ID of the CPU/Core
          * @param nCpu    Number of cores used
          * @param dim     Dimensions
          * @param split   Dimension to split
          */
         SingleSplitting(const int id, const int nCpu, const ArrayI& dim, Splitting::Locations::Id split);

         /**
          * @brief Destructor
          */
         virtual ~SingleSplitting();

         /**
          * @brief Check if factorisation is applicable to scheme
          */
         virtual bool applicable() const;
         
      protected:
         /**
          * @brief Split ith dimension transform
          *
          * @param transId Split the ith dimension
          * @param cpuId   ID of the CPU
          * @param status  Status output
          */
         virtual SharedTransformResolution splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status);

         /**
          * @brief Select the transform grouper
          *
          * \mhdTodo Grouper selection is not correctly implemented
          */
         void selectGrouper();

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         virtual Array computeScore(SharedResolution spResolution);

      private:
         /**
          * @brief Dimension to split
          */
         Splitting::Locations::Id mSplit;
   };

}
}

#endif // QUICC_PARALLEL_SINGLESPLITTING_HPP
