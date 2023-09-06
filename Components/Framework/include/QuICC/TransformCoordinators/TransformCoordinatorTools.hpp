/**
 * @file TransformCoordinatorTools.hpp
 * @brief Implementation of transform coordinator tools
 */

#ifndef QUICC_TRANSFORM_TRANSFORMCOORDINATORTOOLS_HPP
#define QUICC_TRANSFORM_TRANSFORMCOORDINATORTOOLS_HPP

// System includes
//
#include <map>

// Project includes
//
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/TransformGroupers/IForwardGrouper.hpp"
#include "QuICC/TransformGroupers/IBackwardGrouper.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/Timers/StageTimer.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief Implementation of transform coordinator tools
    */
   class TransformCoordinatorTools
   {
      public:
         /**
          * @brief Compute number of packs
          *
          * @param packs            Packs
          * @param spFwdGrouper     Forward transform communication grouper
          * @param spBwdGrouper     Backward transform communication grouper
          * @param forwardTrees     Transform integrator trees
          * @param backwardTrees    Transform projector trees
          * @param keys             Keys for tree maps
          * @param spRes            Shared resolution
          */
         static void computePacks(std::vector<ArrayI>& packs, SharedIForwardGrouper spFwdGrouper, SharedIBackwardGrouper spBwdGrouper, const std::map<int, std::vector<Transform::TransformTree> >& forwardTree, const std::map<int, std::vector<Transform::TransformTree> >& backwardTree, const std::set<int>& keys, SharedResolution spRes);

         /**
          * @brief Initialise the transform coordinator
          *
          * @param rCoord           Transform coordinator
          * @param spFwdGrouper     Forward transform communication grouper
          * @param spBwdGrouper     Backward transform communication grouper
          * @param packs            pack sizes
          * @param spRes            Shared resolution
          * @param runOptions       Available run options map
          */
         static void init(TransformCoordinatorType& rCoord, SharedIForwardGrouper spFwdGrouper, SharedIBackwardGrouper spBwdGrouper, const std::vector<ArrayI>& packs, SharedResolution spRes, const std::map<std::size_t,NonDimensional::SharedINumber>& runOptions);

      protected:

      private:
         /**
          * @brief Add set of packs as ArrayI to packs
          */
         static void addSet(std::vector<ArrayI>& packs, const std::set<int>& packSet);

         /// Tag type for fillSet, stage 1D
         struct stage_1D {};
         /// Tag type for fillSet, stage 2D
         struct stage_2D {};

         /**
          * @brief Fill set with pack size
          */
         template <typename TStage, typename TGrouper> static void fillSet(std::set<int>& packSet, TGrouper spGrouper, const std::vector<Transform::TransformTree>& tree);
         /**
          * @brief Constructor
          */
         TransformCoordinatorTools() = default;

         /**
          * @brief Destructor
          */
         ~TransformCoordinatorTools() = default;
   };

   template <typename TStage, typename TGrouper> void TransformCoordinatorTools::fillSet(std::set<int>& packSet, TGrouper spGrouper, const std::vector<Transform::TransformTree>& tree)
   {
      ArrayI p;
      if constexpr(std::is_same_v<TStage,stage_1D>)
      {
         p = spGrouper->packs1D(tree);
      }
      else if constexpr(std::is_same_v<TStage,stage_2D>)
      {
         p = spGrouper->packs2D(tree);
      }
      else
      {
         static_assert( sizeof(TStage) && false, "unknown stage type");
      }

      for(int i = 0; i < p.size(); i++)
      {
         packSet.insert(p(i));
      }
   }

}
}

#endif // QUICC_TRANSFORM_TRANSFORMCOORDINATORTOOLS_HPP
