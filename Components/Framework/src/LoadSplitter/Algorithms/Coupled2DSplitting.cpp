/** 
 * @file Coupled2DSplitting.cpp
 * @brief Source of the implementation of a load splitting algorithm for coupled 2D matrices
 */

// System includes
//
#include <utility>
#include <queue>
#include <set>
#include <map>

// External includes
//

// Class include
//
#include "QuICC/LoadSplitter/Algorithms/Coupled2DSplitting.hpp"

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/SplittingTools.hpp"

namespace QuICC {

namespace Parallel {

   Coupled2DSplitting::Coupled2DSplitting(const int id, const int nCpu, const ArrayI& dim)
      : SplittingAlgorithm(id, nCpu, dim, Splitting::Algorithms::COUPLED2D)
   {
      // Initialise the NCpu factors
      this->initFactors(1);
   }

   Coupled2DSplitting::~Coupled2DSplitting()
   {
   }

   bool Coupled2DSplitting::applicable() const
   {
      bool status = true;

      // Only works with at least 2 dimensions
      status = (status && (this->dims() > 1));

      // Check for scheme specific conditions
      status = (status && this->mspScheme->applicable());

      return status;
   }

   SharedTransformResolution  Coupled2DSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status)
   {
      // Get size of the splittable dimension(s)
      int tot = this->mspScheme->splittableTotal(transId, Splitting::Locations::COUPLED2D);

      // Build a simple balanced split
      ArrayI ids(1);
      ArrayI n0(1);
      ArrayI nN(1);
      SplittingTools::balancedSplit(n0(0), nN(0), tot, this->factor(0), cpuId);

      // Storage for the forward 1D indexes
      std::vector<ArrayI>  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<ArrayI>  bwd1D;
      // Storage for the 2D indexes
      std::vector<ArrayI>  idx2D;
      // Storage for the 3D indexes
      ArrayI  idx3D;

      // Compute the indexes
      ids(0) = cpuId;
      status = this->mspScheme->fillIndexes(transId, fwd1D, bwd1D, idx2D, idx3D, ids, this->factors(), n0, nN, Splitting::Locations::COUPLED2D);

      // Create TransformResolution object
      auto spTraRes = std::make_shared<TransformResolution>(fwd1D, bwd1D, idx2D, idx3D);
      return spTraRes;
   }

   void Coupled2DSplitting::selectGrouper()
   {
      // SINGLE1D or TRANSFORM grouper setup
      #if defined QUICC_TRANSGROUPER_SINGLE1D
         this->mGrouper = Splitting::Groupers::SINGLE1D;
      #else
         this->mGrouper = Splitting::Groupers::EQUATION;
      #endif //defined QUICC_TRANSGROUPER_SINGLE1D
   }

   Array Coupled2DSplitting::computeScore(SharedResolution spResolution)
   {
      // Initialise the score
      Array details(4);
      details(0) = 100;

      // Multiply by communication score
      ArrayI comm;
      details(1) = this->communicationScore(spResolution, comm);

      // Multiply by load balancing score
      Array balance = this->mspScheme->loadWeights();
      details(2) = this->balancingScore(spResolution, balance);

      // Use additional memory related weighting
      details(3) = this->mspScheme->memoryScore(spResolution);

      // Select best transform grouper algorithm
      this->selectGrouper();

      return details;
   }

}
}
