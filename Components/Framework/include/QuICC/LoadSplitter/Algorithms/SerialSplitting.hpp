/** 
 * @file SerialSplitting.hpp
 * @brief Implementation of a serial "load splitting" algorithm
 */

#ifndef QUICC_PARALLEL_SERIALSPLITTING_HPP
#define QUICC_PARALLEL_SERIALSPLITTING_HPP

// System includes
//
#include <utility>
#include <vector>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "QuICC/Resolutions/TransformResolution.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of a serial "load splitting" algorithm
    */
   class SerialSplitting: public SplittingAlgorithm
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id   ID of the CPU
          * @param nCpu Number of CPUs used
          * @param dim  Dimensions
          */
         SerialSplitting(const int id, const int nCpu, const ArrayI& dim);

         /**
          * @brief Destructor
          */
         ~SerialSplitting() = default;

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
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         Array computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp) final;

      private:
   };

}
}

#endif // QUICC_PARALLEL_SERIALSPLITTING_HPP
