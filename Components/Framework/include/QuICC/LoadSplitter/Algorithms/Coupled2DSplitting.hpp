/**
 * @file Coupled2DSplitting.hpp
 * @brief Implementation of a load splitting algorithm for coupled 2D matrices (only slowest dimension can be split for serial solver)
 */

#ifndef QUICC_PARALLEL_COUPLED2DSPLITTING_HPP
#define QUICC_PARALLEL_COUPLED2DSPLITTING_HPP

// System includes
//
#include <utility>
#include <vector>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "QuICC/Resolutions/TransformResolution.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of a load splitting algorithm for coupled 2D matrices (only slowest dimension can be split for serial solver)
    */
   class Coupled2DSplitting: public SplittingAlgorithm
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id      ID of the CPU/Core
          * @param nCpu    Number of cores used
          * @param dim     Dimensions
          */
         Coupled2DSplitting(const int id, const int nCpu, const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~Coupled2DSplitting();

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
   };

}
}

#endif // QUICC_PARALLEL_COUPLED2DSPLITTING_HPP
