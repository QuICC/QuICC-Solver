/**
 * @file TransformTreeTools.cpp
 * @brief Source of the implementation of tools to work with projector transform trees
 */

// Debug includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// Configuration includes
//

// System includes
//
#include <set>
#include <vector>
#include <map>
#include <tuple>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/TransformConfigurators/TransformTreeTools.hpp"

// Project includes
//
#include "QuICC/Arithmetics/Add.hpp"
#include "QuICC/Arithmetics/Set.hpp"
#include "QuICC/Arithmetics/SetNeg.hpp"
#include "QuICC/Arithmetics/Sub.hpp"
#include "QuICC/Arithmetics/None.hpp"
#include "QuICC/Io/Xml/GxlWriter.hpp"

namespace QuICC {

namespace Transform {

   void TransformTreeTools::generateTrees(std::vector<TransformTree>& rTrees, const std::map<std::size_t, std::vector<TransformPath> >& branches, const TransformDirection::Id dir, const std::string& prepend)
   {
      std::string dirName;
      if(dir == TransformDirection::FORWARD)
      {
         dirName = "forward";
      } else
      {
         dirName = "backward";
      }

      if(!prepend.empty())
      {
         dirName = prepend + "_" + dirName;
      }

      //
      // Output transform branches as a graph
      //
      Io::Xml::GxlWriter gxlPath(dirName + "_transform_paths");
      gxlPath.init();
      gxlPath.graphTransformPath(branches, dir);
      gxlPath.write();
      gxlPath.finalize();

      //
      // First stage: Collapse calculations with same input into same tree branch
      //
      buildTrees(rTrees, branches);

      //
      // Second stage: Set proper recovery and hold flags on the tree nodes
      //
      finalizeTrees(rTrees);

      //
      // Output transfrom tree as a graph
      //
      Io::Xml::GxlWriter gxlTree(dirName + "_transform_trees");
      gxlTree.init();
      gxlTree.graphTransformTree(rTrees, dir);
      gxlTree.write();
      gxlTree.finalize();

      #if defined QUICC_OPTIMIZE_TREE
      //
      // Third stage: Combine arithmetic operations into single tree branch
      //
      optimizeTrees(rTrees);

      //
      // Output transfrom tree as a graph
      //
      Io::Xml::GxlWriter gxlOpti(dirName + "_transform_optimized_trees");
      gxlOpti.init();
      gxlOpti.graphTransformTree(rTrees, dir);
      gxlOpti.write();
      gxlOpti.finalize();
      #endif //defined QUICC_OPTIMIZE_TREE
   }

   void TransformTreeTools::growTree(std::map<std::size_t, std::vector<TransformPath> >::const_iterator nameIt, std::set<int>::const_iterator compIt, const int dim, std::vector<std::size_t>& path, TransformTreeEdge &rPrev)
   {
      // Typedefs to simplify notation
      typedef std::map<std::size_t, int> OpMap;
      typedef std::multimap<std::size_t, std::tuple<std::vector<int>, FieldType::Id, std::size_t> > OpMulti;

      bool endRecursion = false;

      // Extract unique operators
      OpMap op;
      std::pair<OpMap::iterator,bool> opPairIt;
      OpMulti multi;
      for(auto branchIt = nameIt->second.cbegin(); branchIt != nameIt->second.cend(); ++branchIt)
      {
         // Check path
         bool isSame = true;
         for(int i = 0; i < dim; i++)
         {
            isSame = isSame && (branchIt->edge(i).opId() == path.at(i));
         }

         // Count path multiplicity
         if(branchIt->startId() == *compIt && isSame)
         {
            opPairIt = op.insert(std::make_pair(branchIt->edge(dim).opId(), 0));
            opPairIt.first->second += 1;

            // Check for end of recursion
            if(dim == branchIt->nEdges()-1)
            {
               multi.insert(std::make_pair(branchIt->edge(dim).opId(), std::make_tuple(branchIt->edge(dim).outId(), branchIt->fieldId(), branchIt->edge(dim).arithId())));
               endRecursion = true;
            }
         }
      }

      // Check for end of recursion
      if(endRecursion)
      {
         // Create edges
         for(auto opIt = op.cbegin(); opIt != op.cend(); ++opIt)
         {
            std::pair<OpMulti::const_iterator, OpMulti::const_iterator> multiRange = multi.equal_range(opIt->first);
            for(auto multiIt = multiRange.first; multiIt != multiRange.second; ++multiIt)
            {
               TransformTreeEdge &rNext = rPrev.addEdge(opIt->first, opIt->second);

               rNext.setEnd(std::get<0>(multiIt->second), std::get<1>(multiIt->second), std::get<2>(multiIt->second));
            }
         }

      } else
      {
         // Create edges
         for(auto opIt = op.cbegin(); opIt != op.cend(); ++opIt)
         {
            TransformTreeEdge &rNext = rPrev.addEdge(opIt->first, opIt->second);

            // Add last element to path
            path.push_back(opIt->first);

            // Grow the tree
            growTree(nameIt, compIt, dim+1, path, rNext);

            // Remove last element from path
            path.pop_back();
         }
      }
   }

   void TransformTreeTools::buildTrees(std::vector<TransformTree>& rTrees, const std::map<std::size_t, std::vector<TransformPath> >& branches)
   {
      //
      // Construct the trees from the transform paths
      //
      for(auto nameIt = branches.cbegin(); nameIt != branches.cend(); ++nameIt)
      {
         // Generate list of unique components
         std::set<int> components;
         for(auto branchIt = nameIt->second.cbegin(); branchIt != nameIt->second.cend(); ++branchIt)
         {
            components.insert(branchIt->startId());
         }

         // Initialise component trees
         std::vector<std::size_t> path;
         for(auto compIt = components.cbegin(); compIt != components.cend(); ++compIt)
         {
            // Initialize tree
            rTrees.push_back(TransformTree(nameIt->first, *compIt, nameIt->second.at(0).nEdges()));

            // Grow the tree recursively
            path.clear();
            growTree(nameIt, compIt, 0, path, rTrees.back().rRoot());
         }
      }
   }

   void TransformTreeTools::setInputInfoEdge(TransformTreeEdge& edge)
   {
      // Get range of edges
      TransformTreeEdge::EdgeType_range rangeIt = edge.rEdgeRange();

      // Initialize recover and hold
      int recover = 0;
      int hold = std::distance(rangeIt.first, rangeIt.second) - 1;

      // Loop over edges
      for(auto edgeIt = rangeIt.first; edgeIt != rangeIt.second; ++edgeIt, ++recover, --hold)
      {
         // Set recover and hold info
         edgeIt->setInputInfo(recover, hold);

         // Recursively go to next level
         setInputInfoEdge(*edgeIt);
      }
   }

   void TransformTreeTools::setOutputOrder(TransformTreeEdge& edge, TransformTreeTools::SetOutIds& outIds)
   {
      // Get range of edges
      TransformTreeEdge::EdgeType_range rangeIt = edge.rEdgeRange();

      // Check for last recurrence
      if(std::distance(rangeIt.first->edgeRange().first, rangeIt.first->edgeRange().second) > 0)
      {
         // Loop over edges
         for(auto edgeIt = rangeIt.first; edgeIt != rangeIt.second; ++edgeIt)
         {
            // Recursively go to next level
            setOutputOrder(*edgeIt, outIds);
         }

      // Reached last level
      } else
      {
         // Loop over edges
         for(auto edgeIt = rangeIt.first; edgeIt != rangeIt.second; ++edgeIt)
         {
            // Try to add to unique output fields
            std::pair<SetOutIds::iterator,bool> outIt = outIds.insert(std::make_pair(edgeIt->fieldId(), edgeIt->outIds()));

            // If first entry, replace arithmetics with Set or SetNeg
            if(outIt.second)
            {
               if(edgeIt->arithId() == Arithmetics::Add::id())
               {
                  edgeIt->setArithId(Arithmetics::Set::id());

               } else if(edgeIt->arithId() == Arithmetics::Sub::id())
               {
                  edgeIt->setArithId(Arithmetics::SetNeg::id());
               } else
               {
                  throw std::logic_error("Unknown arithmetic operation during transform tree construction!");
               }
            }
         }
      }
   }

   void TransformTreeTools::finalizeTrees(std::vector<TransformTree>& rTrees)
   {
      // Loop over all trees
      SetOutIds outIds;
      std::size_t curName = 0;
      for(auto treeIt = rTrees.begin(); treeIt != rTrees.end(); ++treeIt)
      {
         if(curName != treeIt->name())
         {
            outIds.clear();
            curName = treeIt->name();
         }

         // Recursively set hold and recovery flags
         setInputInfoEdge(treeIt->rRoot());

         // Recursively order operations by replacing first by Set or SetNeg
         setOutputOrder(treeIt->rRoot(), outIds);
      }
   }

   void TransformTreeTools::optimizeOutputCount(TransformTreeEdge& edge, TransformTreeTools::MapOptIds& optIds, SetOutIds& outIds, int& counter)
   {
      // Get range of edges
      TransformTreeEdge::EdgeType_range rangeIt = edge.rEdgeRange();

      // Check for last recurrence
      if(std::distance(rangeIt.first->edgeRange().first, rangeIt.first->edgeRange().second) > 0)
      {
         // Loop over edges
         for(auto edgeIt = rangeIt.first; edgeIt != rangeIt.second; ++edgeIt)
         {
            // Recursively go to next level
            optimizeOutputCount(*edgeIt, optIds, outIds, counter);
         }

      // Reached last level
      } else
      {
         // Loop over edges
         for(auto edgeIt = rangeIt.first; edgeIt != rangeIt.second; ++edgeIt)
         {
            // Add to set of output keys
            OutKey outKey = std::make_pair(edgeIt->fieldId(), edgeIt->outIds());
            std::pair<SetOutIds::iterator,bool> outIt = outIds.insert(outKey);

            // Add to map or increase count
            OptKey optKey = std::make_tuple(edgeIt->opId(), edgeIt->fieldId(), edgeIt->outIds());
            std::pair<MapOptIds::iterator,bool> optIt = optIds.insert(std::make_pair(optKey, std::make_tuple(0,counter,outIt.second)));
            std::get<0>(optIt.first->second) += 1;

            // Increment counter used as index
            counter++;
         }
      }
   }

   void TransformTreeTools::pruneOutput(TransformTreeEdge& edge, TransformTreeTools::MapOptIds& optIds)
   {
      // Get range of edges
      TransformTreeEdge::EdgeType_range rangeIt = edge.rEdgeRange();

      // Check for last recurrence
      if(std::distance(rangeIt.first->edgeRange().first, rangeIt.first->edgeRange().second) > 0)
      {
         // Loop over edges
         for(auto edgeIt = rangeIt.first; edgeIt != rangeIt.second; ++edgeIt)
         {
            // Recursively go to next level
            pruneOutput(*edgeIt, optIds);
         }

      // Reached last level
      } else
      {
         // Algorithm only works for single branch
         bool hasSingleLeaf = (std::distance(edge.edgeRange().first, edge.edgeRange().second) == 1);

         // Get edge iterator and make key
         TransformTreeEdge::EdgeType_iterator edgeIt = edge.rEdgeRange().first;
         OptKey optKey = std::make_tuple(edgeIt->opId(), edgeIt->fieldId(), edgeIt->outIds());

         // Setup first element of combined field
         if(std::get<0>(optIds.find(optKey)->second) > 1)
         {
            // Make counter negative to identify them later on
            std::get<0>(optIds.find(optKey)->second) = -std::get<0>(optIds.find(optKey)->second) + 1;

            // Convert to set
            std::size_t arithId;
            if(edgeIt->arithId() == Arithmetics::Add::id())
            {
               arithId = Arithmetics::Set::id();
            } else if(edgeIt->arithId() == Arithmetics::Sub::id())
            {
               arithId = Arithmetics::SetNeg::id();
            } else
            {
               arithId = edgeIt->arithId();
            }

            edge.setCombinedInfo(arithId, -1, std::get<1>(optIds.find(optKey)->second));
            edge.setArithId(Arithmetics::None::id());
            edgeIt = edge.delEdge(edgeIt);

         // Setup other elements in combination
         } else if(std::get<0>(optIds.find(optKey)->second) < 0)
         {
            if(!hasSingleLeaf)
            {
               throw std::logic_error("Tree optimization algorithm doesn't work on this tree!");
            }

            if(edgeIt->arithId() == Arithmetics::Set::id() || edgeIt->arithId() == Arithmetics::SetNeg::id())
            {
               throw std::logic_error("The tree to optimize is not properly setup!");
            }
            std::size_t arithId = edgeIt->arithId();
            edge.setArithId(Arithmetics::None::id());

            int holdId;
            // Setup intermediate elements
            if(std::get<0>(optIds.find(optKey)->second) < -1)
            {
               holdId = std::get<1>(optIds.find(optKey)->second);
               edgeIt = edge.delEdge(edgeIt);

            // Setup last element
            } else
            {
               holdId = -1;
               if(std::get<2>(optIds.find(optKey)->second))
               {
                  edgeIt->setArithId(Arithmetics::Set::id());
               } else
               {
                  edgeIt->setArithId(Arithmetics::Add::id());
               }
            }
            edge.setCombinedInfo(arithId, std::get<1>(optIds.find(optKey)->second), holdId);
            std::get<0>(optIds.find(optKey)->second) += 1;
         }
      }
   }

   void TransformTreeTools::optimizeTrees(std::vector<TransformTree>& rTrees)
   {
      // Loop over all trees
      if(rTrees.size() > 0)
      {
         SetOutIds outIds;
         MapOptIds optIds;
         std::size_t curName = rTrees.begin()->name();
         std::vector<TransformTree>::iterator subBegin = rTrees.begin();
         int counter = 0;
         for(auto treeIt = rTrees.begin(); treeIt != rTrees.end(); ++treeIt)
         {
            if(curName != treeIt->name())
            {
               curName = treeIt->name();
               if(optIds.size() > 0)
               {
                  for(auto subIt = subBegin; subIt != treeIt; ++subIt)
                  {
                     pruneOutput(subIt->rRoot(), optIds);
                  }
                  subBegin = treeIt;
               }
               outIds.clear();
               optIds.clear();
            }

            // Recursively set hold and recovery flags
            optimizeOutputCount(treeIt->rRoot(), optIds, outIds, counter);
         }

         // Prune last field
         if(optIds.size() > 0)
         {
            for(auto subIt = subBegin; subIt != rTrees.end(); ++subIt)
            {
               pruneOutput(subIt->rRoot(), optIds);
            }
         }
      }
   }
}
}
