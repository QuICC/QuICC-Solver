/** 
 * @file TubularSplitting.cpp
 * @brief Source of the implementation of the double load splitting algorithm, aka "Tubular" splitting
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/LoadSplitter/Algorithms/TubularSplitting.hpp"

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/SplittingTools.hpp"

namespace QuICC {

namespace Parallel {

   TubularSplitting::TubularSplitting(const int id, const int nCpu, const ArrayI& dim)
      : SplittingAlgorithm(id, nCpu, dim, Splitting::Algorithms::TUBULAR)
   {
      // Initialise the NCpu factors
      this->initFactors(2);
   }

   TubularSplitting::~TubularSplitting()
   {
   }

   bool TubularSplitting::applicable() const
   {
      bool status = true;

      // Check that all three dimensions are splittable by factors
      status = (status && (this->mspScheme->splittableTotal(Dimensions::Transform::TRA1D, Splitting::Locations::SECOND) >= this->factor(1)));
      status = (status && (this->mspScheme->splittableTotal(Dimensions::Transform::TRA1D, Splitting::Locations::BOTH) >= this->factor(0)));
      status = (status && (this->mspScheme->splittableTotal(Dimensions::Transform::TRA2D, Splitting::Locations::BOTH) >= this->factor(1)));
      status = (status && (this->mspScheme->splittableTotal(Dimensions::Transform::TRA2D, Splitting::Locations::FIRST) >= this->factor(0)));
      status = (status && (this->mspScheme->splittableTotal(Dimensions::Transform::TRA3D, Splitting::Locations::FIRST) >= this->factor(0)));
      status = (status && (this->mspScheme->splittableTotal(Dimensions::Transform::TRA3D, Splitting::Locations::BOTH) >= this->factor(1)));

      return status;
   }

   SharedTransformResolution  TubularSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status)
   {
      // Storage for the forward 1D indexes
      std::vector<ArrayI>  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<ArrayI>  bwd1D;
      // Storage for the 2D indexes
      std::vector<ArrayI>  idx2D;
      // Storage for the 3D indexes
      ArrayI  idx3D;

      // Create arrays for the IDs
      ArrayI ids(2);
      ids(0) = SplittingTools::groupId(this->factors(), 0, cpuId);
      ids(1) = SplittingTools::groupId(this->factors(), 1, cpuId);

      // Create start indexes and length
      ArrayI n0;
      ArrayI nN;

      if(transId == Dimensions::Transform::TRA1D)
      {
         // Get size of the dimension that needs to be compatible with previous step
         int tot = this->mspScheme->splittableTotal(transId, Splitting::Locations::SECOND);

         // Build a balanced split (has to be the same as in TRA2D)
         int c0, cN;
         SplittingTools::balancedSplit(c0, cN, tot, this->factor(1), ids(1));

         // Get the remaining splittable size
         tot = cN*this->mspScheme->splittableTotal(transId, Splitting::Locations::BOTH);

         // Build a simple balanced split over remaining modes
         int r0, rN;
         SplittingTools::balancedSplit(r0, rN, tot, this->factor(0), ids(0));

         // 
         // The compatible dimension is the second dimension, we will need to reorder the indexes
         //

         // Get bottom row offset
         int b0 = r0 % cN;

         // Get the top row length
         int tN = (rN-(cN-b0)) % cN;
         // WARNING:  tN = 0 means that it is a full row!
         if(tN == 0)
         {
            tN = cN;
         }

         // Compute shifted offset
         int s0 = r0/cN;
         // ... and shifted size (+1 because of bottom row)
         int sN = static_cast<int>(std::ceil(static_cast<double>(rN-(cN-b0))/static_cast<double>(cN))) + 1;

         // Create start indexes and length
         n0.resize(sN+1);
         nN.resize(sN+1);
         n0.setConstant(0);
         nN.setConstant(0);
         n0(0) = s0;
         nN(0) = sN;

         // General setup
         n0.tail(sN).setConstant(c0);
         nN.tail(sN).setConstant(cN);

         // Special treatment for bottom row
         n0(1) = (c0 + b0);
         nN(1) = cN - b0;

         // Special treatment for top row
         n0.tail(1).setConstant(c0);
         nN.tail(1).setConstant(tN);

      } else if(transId == Dimensions::Transform::TRA2D)
      {
         // Create start indexes and length
         n0.resize(2);
         nN.resize(2);

         // Get size of the splittable dimension(s)
         int tot = this->mspScheme->splittableTotal(transId, Splitting::Locations::BOTH);

         // Build a simple balanced split
         SplittingTools::balancedSplit(n0(0), nN(0), tot, this->factor(1), ids(1));

         // Get size of the splittable dimension(s)
         tot = this->mspScheme->splittableTotal(transId, Splitting::Locations::FIRST);

         // Compute a balanced splitting
         SplittingTools::balancedSplit(n0(1), nN(1), tot, this->factor(0), ids(0));

      } else if(transId == Dimensions::Transform::TRA3D)
      {
         // Get size of the splittable dimension(s)
         int tot = this->mspScheme->splittableTotal(transId, Splitting::Locations::FIRST);

         // Compute a balanced splitting
         int t0, tN;
         SplittingTools::balancedSplit(t0, tN, tot, this->factor(0), ids(0));

         // Create start indexes and length
         n0.resize(tN+1);
         nN.resize(tN+1);
         n0.setConstant(0);
         nN.setConstant(0);
         n0(0) = t0;
         nN(0) = tN;

         // Get total for second dimension
         tot = nN(0)*this->mspScheme->splittableTotal(transId, Splitting::Locations::BOTH);

         // Get balanced splitting for second direction
         SplittingTools::balancedSplit(t0, tN, tot, this->factor(1), ids(1));

         // Compute starting point of grid points
         for(int i = 0; i < t0; ++i)
         {
            n0((i % nN(0)) + 1) += 1;
         }

         // Get smallest th0 to use as offset
         t0 = 0;
         int n0Min = n0(1);
         for(int r = 1; r < nN(0); ++r)
         {
            if(n0Min > n0(r+1))
            {
               t0 = r;
               n0Min = n0(r+1);
            }
         }

         // Compute optimal number of grid points
         for(int i = 0; i < tN; ++i)
         {
            nN(((i + t0) % nN(0)) + 1) += 1;
         }
      }

      // Compute the indexes
      status = this->mspScheme->fillIndexes(transId, fwd1D, bwd1D, idx2D, idx3D, ids, this->factors(), n0, nN, Splitting::Locations::BOTH);

      // Create TransformResolution object
      auto spTraRes = std::make_shared<TransformResolution>(fwd1D, bwd1D, idx2D, idx3D);
      return spTraRes;
   }

   void TubularSplitting::selectGrouper(const Splitting::Groupers::Id selected)
   {
      this->mGrouper = selected;
   }

   Array TubularSplitting::computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp)
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
