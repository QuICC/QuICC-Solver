/** 
 * @file SerialSplitting.cpp
 * @brief Source of the implementation of a serial code "load splitting"
 */

// System includes
//

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/SerialSplitting.hpp"

namespace QuICC {

namespace Parallel {

   SerialSplitting::SerialSplitting(const int id, const int nCpu, const ArrayI& dim)
      : SplittingAlgorithm(id, nCpu, dim, Splitting::Algorithms::SERIAL)
   {
      // Initialise the NCpu factors
      this->initFactors(1);
   }

   bool SerialSplitting::applicable() const
   {
      // As long as there is a single CPU it should be applicable
      bool status = this->nCpu() == 1;

      return status;
   }

   SharedTransformResolution  SerialSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status)
   {
      // Create arrays for the IDs and bins
      std::vector<int> ids = {0,0};
      std::vector<int> bins = {1,1};

      // Storage for the forward 1D indexes
      std::vector<std::vector<std::vector<int> > >  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<std::vector<std::vector<int> > >  bwd1D;
      // Storage for the 2D indexes
      std::vector<std::vector<int> >  idx2D;
      // Storage for the 3D indexes
      std::vector<int>  idx3D;

      // Compute the indexes
      status = this->mspScheme->fillIndexes(transId, fwd1D, bwd1D, idx2D, idx3D, ids, bins);

      // Create TransformResolution object
      auto spTraRes = std::make_shared<TransformResolution>(fwd1D, bwd1D, idx2D, idx3D);
      return spTraRes;
   }

   Array SerialSplitting::computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp)
   {
      return Array::Constant(4,1.0);
   }
}
}
