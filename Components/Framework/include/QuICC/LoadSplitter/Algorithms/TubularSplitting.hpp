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

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
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
          * @param split   Dimension to split
          */
         TubularSplitting(const int id, const int nCpu, const ArrayI& dim, Splitting::Locations::Id split);

         /**
          * @brief Destructor
          */
         virtual ~TubularSplitting();

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
          */
         void selectGrouper(const Splitting::Groupers::Id selected);

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         virtual Array computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp);

      private:
         /**
          * @brief Compute a balanced split for 2D data distribution
          */
         void balanced2DSplit(ArrayI& n0, ArrayI& nN, const Dimensions::Transform::Id transId, const ArrayI& bins, const ArrayI& ids, const std::vector<Splitting::Locations::Id>& locs, const bool combine, const bool allowEmpty = false);

         /**
          * @brief Dimension to split
          */
         Splitting::Locations::Id mSplit;
   };

}
}

#endif // QUICC_PARALLEL_TUBULARSPLITTING_HPP
