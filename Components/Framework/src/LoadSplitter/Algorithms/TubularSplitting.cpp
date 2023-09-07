/**
 * @file TubularSplitting.cpp
 * @brief Source of the implementation of the double load splitting algorithm, aka "Tubular" splitting
 */

// System includes
//

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/TubularSplitting.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingTools.hpp"

namespace QuICC {

namespace Parallel {

   TubularSplitting::TubularSplitting(const int id, const int nCpu, const ArrayI& dim, Splitting::Algorithms::Id algorithm, const std::list<int>& factors)
      : SplittingAlgorithm(id, nCpu, dim, algorithm)
   {
      // Initialise the NCpu factors
      this->initFactors(2);

      if(algorithm == Splitting::Algorithms::SINGLE1D)
      {
         std::list<int> fList = {nCpu, 1};
         this->useFactorization(fList);
      }
      else if(algorithm == Splitting::Algorithms::SINGLE2D)
      {
         std::list<int> fList = {1, nCpu};
         this->useFactorization(fList);
      }
      // Use imposed factorization factors
      else
      {
         this->useFactorization(factors);
      }
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
      // Create arrays for the IDs
      std::vector<int> ids;
      ids.push_back(SplittingTools::groupId(this->factors(), 0, cpuId));
      ids.push_back(SplittingTools::groupId(this->factors(), 1, cpuId));

      // Create bins
      std::vector<int> bins;
      bins.push_back(this->factors()(0));
      bins.push_back(this->factors()(1));

      // Use 1xNcpu splitting for SPECTRAL
      if(!this->mspScheme->sameSpectralOrdering() && transId == Dimensions::Transform::SPECTRAL)
      {
         bins.at(0) = 1;
         bins.at(1) = this->factor(0)*this->factor(1);
         ids.at(0) = 0;
         ArrayI totBins(1);
         totBins << this->factor(0)*this->factor(1);
         ids.at(1) = SplittingTools::groupId(totBins, 0, cpuId);
      }

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

   void TubularSplitting::selectGrouper(const Splitting::Groupers::Id selected)
   {
      // Only split in first transpose
      if(this->factors()(1) == 1 && selected != Splitting::Groupers::SINGLE1D)
      {
         this->mGrouper = Splitting::Groupers::EQUATION;
      }
      // Only split in second transpose
      else if(this->factors()(0) == 1 && selected != Splitting::Groupers::SINGLE2D)
      {
         this->mGrouper = Splitting::Groupers::EQUATION;
      }
      else
      {
         this->mGrouper = selected;
      }
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
