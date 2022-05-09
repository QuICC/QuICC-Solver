/** 
 * @file SerialSplitting.cpp
 * @brief Source of the implementation of a serial code "load splitting"
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/LoadSplitter/Algorithms/SerialSplitting.hpp"

// Project includes
//

namespace QuICC {

namespace Parallel {

   SerialSplitting::SerialSplitting(const int id, const int nCpu, const ArrayI& dim)
      : SplittingAlgorithm(id, nCpu, dim, Splitting::Algorithms::SERIAL)
   {
      // Initialise the NCpu factors
      this->initFactors(1);
   }

   SerialSplitting::~SerialSplitting()
   {
   }

   bool SerialSplitting::applicable() const
   {
      bool status = true;

      // As long as there is a single CPU it should be applicable
      status = this->nCpu() == 1;

      // Check for scheme specific conditions
      status = (status && this->mspScheme->applicable());

      return status;
   }

   SharedTransformResolution  SerialSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status)
   {
      // Storage for the forward 1D indexes
      std::vector<ArrayI>  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<ArrayI>  bwd1D;
      // Storage for the 2D indexes
      std::vector<ArrayI>  idx2D;
      // Storage for the 3D indexes
      ArrayI  idx3D;

      // Compute the indexes
      status = this->mspScheme->fillIndexes(transId, fwd1D, bwd1D, idx2D, idx3D);

      // Create TransformResolution object
      auto spTraRes = std::make_shared<TransformResolution>(fwd1D, bwd1D, idx2D, idx3D);
      return spTraRes;
   }

   Array SerialSplitting::computeScore(SharedResolution spResolution)
   {
      return Array::Constant(4,1.0);
   }
}
}
