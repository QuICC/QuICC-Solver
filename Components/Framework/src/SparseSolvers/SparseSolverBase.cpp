/**
 * @file SparseSolverBase.cpp
 * @brief Implementation of the base for linear solver structures
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/SparseSolvers/SparseSolverBase.hpp"

// Project includes
//

namespace QuICC {

namespace Solver {

   SparseSolverBase::SparseSolverBase(const int start, const std::size_t timeId)
      : mZeroIdx(start), mSolveTiming(timeId), mIsInitialized(false)
   {
      // Safety assert
      assert(start >= 0);
   }

   SparseSolverBase::~SparseSolverBase()
   {
   }

   void SparseSolverBase::addInformation(const SpectralFieldId& id, const int idx, const ArrayI& startRow)
   {
      this->mFieldIds.push_back(id);

      this->mInformation.insert(std::make_pair(id, std::make_pair(idx,startRow)));
   }

   void SparseSolverBase::initStartRow()
   {
      // Iterate over all stored fields to build start matrix
      MatrixI start = MatrixI::Zero(this->mInformation.begin()->second.second.rows(), this->mInformation.size());
      int i = 0;
      for(auto jIt = this->mInformation.begin(); jIt != this->mInformation.end(); ++jIt)
      {
         for(auto iIt = this->mInformation.begin(); iIt != this->mInformation.end(); ++iIt)
         {
            if(jIt->second.first > iIt->second.first)
            {
               start.col(i) += iIt->second.second;
            }
         }

         ++i;
      }

      // Store start rows
      i = 0;
      for(auto jIt = this->mInformation.begin(); jIt != this->mInformation.end(); ++jIt)
      {
         jIt->second.second = start.col(i);
         ++i;
      }
   }

   int SparseSolverBase::startRow(const SpectralFieldId& id, const int i) const
   {
      return this->mInformation.find(id)->second.second(i);
   }

   SparseSolverBase::FieldId_range SparseSolverBase::fieldRange() const
   {
      return std::make_pair(this->mFieldIds.begin(), this->mFieldIds.end());
   }

   bool SparseSolverBase::isInitialized() const
   {
      return this->mIsInitialized;
   }

   void SparseSolverBase::setInitialized()
   {
      this->mIsInitialized = true;
   }

   std::size_t SparseSolverBase::solveTiming() const
   {
      return this->mSolveTiming;
   }
}
}
