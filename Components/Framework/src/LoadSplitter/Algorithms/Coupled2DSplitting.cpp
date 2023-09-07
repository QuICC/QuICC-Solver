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

      return status;
   }

   SharedTransformResolution  Coupled2DSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status)
   {
      // Build a simple balanced split
      std::vector<int > ids = {0,0};
      std::vector<int > bins = {1,1};

      // Storage for the forward 1D indexes
      std::vector<std::vector<std::vector<int> > >  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<std::vector<std::vector<int> > >  bwd1D;
      // Storage for the 2D indexes
      std::vector<std::vector<int> >  idx2D;
      // Storage for the 3D indexes
      std::vector<int>  idx3D;

      // Compute the indexes
      ids.at(0) = cpuId;
      bins.at(0) = this->factors()(0);
      status = this->mspScheme->fillIndexes(transId, fwd1D, bwd1D, idx2D, idx3D, ids, bins);

      // Create TransformResolution object
      auto spTraRes = std::make_shared<TransformResolution>(fwd1D, bwd1D, idx2D, idx3D);
      return spTraRes;
   }

   void Coupled2DSplitting::selectGrouper(const Splitting::Groupers::Id selected)
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

   Array Coupled2DSplitting::computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp)
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
