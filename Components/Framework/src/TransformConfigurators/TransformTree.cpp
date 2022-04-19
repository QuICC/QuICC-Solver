/**
 * @file TransformTree.cpp
 * @brief Source of the implementation of the tranform tree
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/TransformConfigurators/TransformTree.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

   TransformTree::TransformTree(const std::size_t name, const int comp, const int depth)
      : mIsActive(true), mName(name), mComp(comp), mDepth(depth), mRoot(-1, 1)
   {
   }

   TransformTree::~TransformTree()
   {
   }

   bool TransformTree::isActive() const
   {
      return this->mIsActive;
   }

   std::size_t TransformTree::name() const
   {
      return this->mName;
   }

   const TransformTreeEdge& TransformTree::root() const
   {
      return this->mRoot;
   }

   TransformTreeEdge& TransformTree::rRoot()
   {
      return this->mRoot;
   }

   int TransformTree::depth() const
   {
      return this->mDepth;
   }

   int TransformTree::nEdges(const int depth) const
   {
      return this->mRoot.nEdges(depth);
   }

   void TransformTree::disable()
   {
      this->mIsActive = false;
   }

}
}
