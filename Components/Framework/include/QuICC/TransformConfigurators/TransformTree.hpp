/**
 * @file TransformTree.hpp
 * @brief This template describes the complete projection tree for 2D space
 */

#ifndef QUICC_TRANSFORM_TRANSFORMTREE_HPP
#define QUICC_TRANSFORM_TRANSFORMTREE_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/TransformConfigurators/TransformTreeEdge.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This template describes the complete projection tree for 2D space
    */
   class TransformTree
   {
      public:
         /**
          * @brief Contructor for operation
          */
         TransformTree(const std::size_t name, const int comp, const int depth);

         /**
          * @brief Destructor
          */
         ~TransformTree();

         /**
          * @brief Is tree active?
          */
         bool isActive() const;

         /**
          * @brief Get physical name
          */
         std::size_t name() const;

         /**
          * @brief Get field component
          */
         template <typename TId> TId comp() const;

         /**
          * @brief depth of tree = number of transforms
          */
         int depth() const;

         /**
          * @brief Number of edges at given depth
          */
         int nEdges(const int depth) const;

         /**
          * @brief Get root of the tree
          */
         const TransformTreeEdge& root() const;

         /**
          * @brief Enable tree
          */
         void enable();

         /**
          * @brief Disable tree
          */
         void disable();

         /**
          * @brief Set root of the tree
          */
         TransformTreeEdge& rRoot();

      protected:

      private:
         /**
          * @brief Is tree active?
          */
         bool mIsActive;

         /**
          * @brief Physical field name
          */
         std::size_t mName;

         /**
          * @brief Root field component
          */
         int mComp;

         /**
          * @brief Depth of tree
          */
         int mDepth;

         /**
          * @brief Root of the tree
          */
         TransformTreeEdge mRoot;
   };

   template <typename TId> inline TId TransformTree::comp() const
   {
      return static_cast<TId>(this->mComp);
   }

}
}

#endif // QUICC_TRANSFORM_TRANSFORMTREE_HPP
