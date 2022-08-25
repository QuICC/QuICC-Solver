/** 
 * @file SingleSplitting.cpp
 * @brief Source of the implementation of a single load splitting algorithm
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
#include "QuICC/LoadSplitter/Algorithms/SingleSplitting.hpp"

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/SplittingTools.hpp"

namespace QuICC {

namespace Parallel {

   SingleSplitting::SingleSplitting(const int id, const int nCpu, const ArrayI& dim, Splitting::Locations::Id split)
      : SplittingAlgorithm(id, nCpu, dim, Splitting::Algorithms::SINGLE1D), mSplit(split)
   {
      // Make sure the right algorithm is set
      if(split == Splitting::Locations::SECOND)
      {
         this->mAlgo = Splitting::Algorithms::SINGLE2D;
      }

      // Initialise the NCpu factors
      this->initFactors(1);
   }

   SingleSplitting::~SingleSplitting()
   {
   }

   bool SingleSplitting::applicable() const
   {
      bool status = true;

      // Only works with at least 2 dimensions
      status = (status && (this->dims() > 1));

      // Make sure split and dimensions are compatible
      if(this->dims() > 2 || this->mSplit == Splitting::Locations::FIRST)
      {
         // All CPUs should have something to do
         int tot;
         for(int i = 0; i < this->dims(); i++)
         {
            tot = this->mspScheme->splittableTotal(static_cast<Dimensions::Transform::Id>(i), this->mSplit);

            status = (status && (tot >= this->factor(0)));
         }
      } else
      {
         status = false;
      }

      // Check for scheme specific conditions
      status = (status && this->mspScheme->applicable());

      return status;
   }

   SharedTransformResolution  SingleSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status)
   {
      // Get size of the splittable dimension(s)
      int tot = this->mspScheme->splittableTotal(transId, this->mSplit);

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
      status = this->mspScheme->fillIndexes(transId, fwd1D, bwd1D, idx2D, idx3D, ids, this->factors(), n0, nN, this->mSplit);

      // Create TransformResolution object
      auto spTraRes = std::make_shared<TransformResolution>(fwd1D, bwd1D, idx2D, idx3D);
      return spTraRes;
   }

   void SingleSplitting::selectGrouper(const Splitting::Groupers::Id selected)
   {
      // Different splitting directions require a different treatment
      if(this->mSplit == Splitting::Locations::FIRST)
      {
         if(selected != Splitting::Groupers::SINGLE1D)
         {
            this->mGrouper = Splitting::Groupers::EQUATION;
         }
         else
         {
            this->mGrouper = selected;
         }
      }
      else
      {
         if(selected != Splitting::Groupers::SINGLE2D)
         {
            this->mGrouper = Splitting::Groupers::EQUATION;
         }
         else
         {
            this->mGrouper = selected;
         }
      }
   }

   Array SingleSplitting::computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp)
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
      this->selectGrouper(grp);

      return details;
   }

}
}
