/** 
 * @file IBuilder.cpp
 * @brief Source of the base for the scheme builder implementations
 */

// Configuration includes
//

// System includes
//
#include <set>
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

   bool IBuilder::sameSpectralOrdering() const
   {
      return true;
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
      for(int i = 0; i < this->dims()+1; i++)
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
      assert(static_cast<int>(transId) < this->mDimensions.size());

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
      assert(static_cast<int>(transId) < this->mDims || transId == Dimensions::Transform::SPECTRAL);

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
         for(const auto& st: descr.structure)
         {
            // Extract the communication group from structure
            std::set<int> filter;
            filter.insert(QuICCEnv().id());
            auto sze = filter.size();

            // Loop to check for nodes only reachable via another node
            for(auto lvl = 0; lvl < 10; lvl++)
            {
               std::set<int> prev = filter;
               for(auto fi: prev)
               {
                  // Get oneway communication nodes
                  for(auto&& p: st.second)
                  {
                     if(p.second == fi)
                     {
                        filter.insert(p.first);
                     }
                  }

                  // Get nodes with direct communication
                  auto seq = st.second.equal_range(fi);
                  for(auto it = seq.first; it != seq.second; ++it)
                  {
                     filter.insert(it->second);
                  }
               }

               if(sze == filter.size())
               {
                  break;
               }
               else
               {
                  sze = filter.size();
               }
            }

            // Convert set to array of CPUs in group
            std::vector<int> groupCpu;
            groupCpu.reserve(filter.size());
            for(auto it = filter.begin(); it != filter.end(); ++it)
            {
               groupCpu.push_back(*it);
            }

            QuICCEnv().addCommunicator(st.first, groupCpu);

            // Synchronize
            QuICCEnv().synchronize();
         }

      #endif //defined QUICC_MPI
   }

}
}
