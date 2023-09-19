/**
 * @file ICosts.hpp
 * @brief Base implementation of a cost based spatial scheme
 */

#ifndef QUICC_SPATIALSCHEME_ICOSTS_HPP
#define QUICC_SPATIALSCHEME_ICOSTS_HPP

// System includes
//
#include <vector>
#include <deque>
#include <queue>
#include <map>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Resolutions/Resolution.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Base implementation of a cost based spatial scheme
    */
   class ICosts
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dims Dimension of the domain
          */
         ICosts(const int dims);

         /**
          * @brief Destructor
          */
         virtual ~ICosts() = default;

         /**
          * @brief Get load balancing weights
          */
         virtual Array loadWeights();

         /**
          * @brief Get memory related score weight
          */
         virtual double memoryScore(SharedResolution spRes);

      protected:
         /**
          * @brief Set transform costs
          */
         virtual void setCosts() = 0;

         /**
          * @brief Set transform scalings
          */
         virtual void setScalings() = 0;

         /**
          * @brief Set transform memory footprint
          */
         virtual void setMemoryScore() = 0;

         /**
          * @brief Set general cost for a transform
          */
         void setCost(MHDFloat c, Dimensions::Transform::Id id);

         /**
          * @brief Set scaling cost for a transform
          */
         void setScaling(MHDFloat c, Dimensions::Transform::Id id);

         /**
          * @brief Set memory cost for a transform
          */
         void setMemory(MHDFloat c, Dimensions::Transform::Id id);

         /**
          * @brief Reset all the load related variables
          */
         void resetLoad();

         /**
          * @brief Update the load sum for the regular load
          */
         void updateLoad(const int parts);

         /**
          * @brief Storage for the list of loads
          */
         std::deque<int>   mLoadList;

         /**
          * @brief Storage for the list of optimal load
          */
         std::queue<int>   mOptimalLoad;

         /**
          * @brief Storage for summed load
          */
         std::vector<int>   mLoad;

         /**
          * @brief Storage for the regularized load
          */
         std::multimap<int, int>   mRegularLoad;

      private:
         /**
          * @brief Cost of the transforms
          *
          * A value smaller than 1 increases the negative effect of inbalance.
          * A value bigger than 1 reduces the negative effect of inbalance.
          */
         Array   mCosts;

         /**
          * @brief Scaling of the transforms
          *
          * A value smaller than 1 increases the negative effect of inbalance.
          * A value bigger than 1 reduces the negative effect of inbalance.
          */
         Array   mScalings;

         /**
          * @brief Memory footprint of the transforms
          *
          * A value smaller than 1 increases the negative effect of inbalance.
          * A value bigger than 1 reduces the negative effect of inbalance.
          */
         Array   mMemory;

   };
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_ICOSTS_HPP
