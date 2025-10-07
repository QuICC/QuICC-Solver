/**
 * @file KaHIPSplitting.hpp
 * @brief Implementation of load distribution using KaHIP library
 */

#ifndef QUICC_PARALLEL_KAHIPSLITTING_HPP
#define QUICC_PARALLEL_KAHIPSLITTING_HPP

// System includes
//
#include <utility>
#include <vector>

// Project includes
//
#include "QuICC/Enums/Splitting.hpp"
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingAlgorithm.hpp"
#include "QuICC/Resolutions/TransformResolution.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of a load distribution algorithm using KaHIP library
    */
   class KaHIPSplitting: public SplittingAlgorithm
   {
      public:
         /**
          * @brief Constructor
          *
          * @param id ID of the CPU/Core
          * @param nCpu Number of cores used
          * @param dim  Dimensions
          * @param algorithm  Splitting algorithm
          * @param factors Imposed CPU factorizations
          */
         KaHIPSplitting(const int id, const int nCpu, const ArrayI& dim, Splitting::Algorithms::Id algorithm, const std::list<int>& factors);

         /**
          * @brief Destructor
          */
         ~KaHIPSplitting() = default;

         /**
          * @brief Check if factorisation is applicable to scheme
          */
         bool applicable() const final;

      protected:
         /**
          * @brief Split ith dimension transform
          *
          * @param transId Split the ith dimension
          * @param cpuId   ID of the CPU
          * @param status  Status output
          */
         SharedTransformResolution splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status) final;

         /**
          * @brief Select the transform grouper
          */
         void selectGrouper(const Splitting::Groupers::Id selected);

         /**
          * @brief Compute the score of the Resolution
          *
          * @param spResolution Shared resolution object
          */
         Array computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp) final;

      private:
         void computeEdgeCut();
         void mapNodes(std::map<std::pair<int,int>,int>& nodes, std::map<std::pair<int,int>,std::vector<int>>* pmodes, std::map<std::pair<int,int>,std::vector<int>>* pgrid, const Dimensions::Transform::Id transId, const int start);

         std::vector<int> mPartition;
         std::map<std::pair<int,int>,int> mMap1Dnodes;
         std::map<std::pair<int,int>,int> mMap2Dnodes;
         std::map<std::pair<int,int>,int> mMap3Dnodes;
   };

}
}

#endif // QUICC_PARALLEL_KAHIPSLITTING_HPP
