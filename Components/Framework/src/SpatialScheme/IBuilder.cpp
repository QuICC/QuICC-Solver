/** 
 * @file IBuilder.cpp
 * @brief Source of the base for the scheme builder implementations
 */

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/IBuilder.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Resolutions/Tools/RegularIndexCounter.hpp"
#include "QuICC/Framework/MpiFramework.hpp"

namespace QuICC {

namespace SpatialScheme {

   void IBuilder::interpretConfigDimensions(ArrayI& rDim)
   {
   }

   void IBuilder::tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      this->tuneMpiResolution(descr);
   }

   IBuilder::IBuilder(const int dims,  const GridPurpose::Id purpose)
      : ICosts(dims), mPurpose(purpose), mDims(dims)
   {
   }

   IBuilder::~IBuilder()
   {
   }

   void IBuilder::addIndexCounter(SharedResolution spRes)
   {
      auto spCounter = std::make_shared<RegularIndexCounter>(spRes->spSim(), spRes->spCpu());

      spRes->setIndexCounter(spCounter);

   }

   GridPurpose::Id IBuilder::purpose() const
   {
      return this->mPurpose;
   }

   void IBuilder::init()
   {
      // Initialise Storage for the dimensions
      for(int i = 0; i < this->dims(); i++)
      {
         this->mDimensions.push_back(ArrayI(this->dims()+1));
      }

      // Initialise the domain's dimensions
      this->setDimensions(); 

      // Set the transform costs
      this->setCosts();

      // Set the transform scalings
      this->setScalings();

      // Set the memory costs
      this->setMemoryScore();
   }

   int IBuilder::dim(const Dimensions::Transform::Id transId, const Dimensions::Data::Id dataId) const
   {
      // Assert for domain dimensions
      assert(static_cast<int>(transId) < this->mDims);

      // Assert for dimension size
      assert(this->mDimensions.at(static_cast<int>(transId)).size() > static_cast<int>(dataId));

      return this->mDimensions.at(static_cast<int>(transId))(static_cast<int>(dataId));
   }

   void IBuilder::setTransformSpace(const ArrayI& dim)
   {
      this->mTransformSpace = dim;
   }

   const ArrayI& IBuilder::getTransformSpace() const
   {
      assert(this->mTransformSpace.size() == this->dims());
      assert(this->mTransformSpace.array().abs().minCoeff() > 0);

      return this->mTransformSpace;
   }

   void IBuilder::setDimension(int n, const Dimensions::Transform::Id transId, const Dimensions::Data::Id dataId)
   {
      // Assert for positive size
      assert(n > 0);

      // Assert for domain dimensions
      assert(static_cast<int>(transId) < this->mDims);

      // Assert for dimension size
      assert(this->mDimensions.at(static_cast<int>(transId)).size() > static_cast<int>(dataId));

      this->mDimensions.at(static_cast<int>(transId))(static_cast<int>(dataId)) = n;
   }

   int IBuilder::dims() const
   {
      return this->mDims;
   }

   void IBuilder::tuneMpiResolution(const Parallel::SplittingDescription& descr)
   {
      #if defined QUICC_MPI
         MpiFramework::initTransformComm(descr.structure.size());
         for(auto vIt = descr.structure.cbegin(); vIt != descr.structure.cend(); ++vIt)
         {
            // Extract the communication group from structure
            std::set<int> filter;
            filter.insert(QuICCEnv().id());
            for(auto it = vIt->equal_range(QuICCEnv().id()).first; it != vIt->equal_range(QuICCEnv().id()).second; ++it)
            {
               filter.insert(it->second);
            }

            // Convert set to array of CPUs in group
            ArrayI groupCpu(filter.size());
            int i = 0;
            for(auto it = filter.begin(); it != filter.end(); ++it)
            {
               groupCpu(i) = *it;
               ++i;
            }

            MpiFramework::addTransformComm(groupCpu);

            // Synchronize
            QuICCEnv().synchronize();
         }

      #endif //defined QUICC_MPI
   }

}
}
