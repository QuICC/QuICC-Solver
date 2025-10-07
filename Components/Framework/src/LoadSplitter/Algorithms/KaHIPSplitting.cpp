/**
 * @file KaHIPSplitting.cpp
 * @brief Source of the implementation of a load splitting algorithm using KaHIP
 */

// System includes
//
#include "kaHIP_interface.h"
#include <iostream>
#include <fstream>
#include <set>

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/KaHIPSplitting.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingTools.hpp"

namespace QuICC {

namespace Parallel {

   namespace details {
   void writeDot(std::string filename, const std::vector<int>& xnodes, const std::vector<int>& xadj, const std::vector<int>& adjncy, const std::vector<int>& part, const int stage = -1);
   void writePartition(const std::string filename, const std::vector<int>& part);
   void writeMetis(const std::string filename, const std::vector<int>& xadj, const std::vector<int>& adjncy, const std::vector<int>& vwgt, const std::vector<int>& adjcwgt);
   }

   KaHIPSplitting::KaHIPSplitting(const int id, const int nCpu, const ArrayI& dim, Splitting::Algorithms::Id algorithm, const std::list<int>& factors)
      : SplittingAlgorithm(id, nCpu, dim, algorithm)
   {
      // Initialise the NCpu factors
      this->initFactors(1);
   }

   bool KaHIPSplitting::applicable() const
   {
      bool status = true;

      // Check that all three dimensions are splittable by factors

      return status;
   }
   
   void KaHIPSplitting::mapNodes(std::map<std::pair<int,int>,int>& nodes, std::map<std::pair<int,int>,std::vector<int>>* pmodes, std::map<std::pair<int,int>,std::vector<int>>* pgrid, const Dimensions::Transform::Id transId, const int start)
   {
      // Create arrays for the IDs and bins
      std::vector<int> ids = {0,0};
      std::vector<int> bins = {1,1};

      // Storage for the forward 1D indexes
      std::vector<std::vector<std::vector<int> > >  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<std::vector<std::vector<int> > >  bwd1D;
      // Storage for the 2D indexes
      std::vector<std::vector<int> >  idx2D;
      // Storage for the 3D indexes
      std::vector<int>  idx3D;

      int status = this->mspScheme->fillIndexes(transId, fwd1D, bwd1D, idx2D, idx3D, ids, bins);
      int nodeId = start;
      for(int j = 0; j < idx3D.size(); j++)
      {
         auto k0 = idx3D.at(j);
         for(int i = 0; i < idx2D.at(j).size(); i++)
         {
            auto k1 = idx2D.at(j).at(i);
            nodes.try_emplace({k0,k1}, nodeId);
            if(pmodes)
            {
               pmodes->try_emplace({k0,k1}, bwd1D.at(j).at(i));
            }
            if(pgrid)
            {
               pgrid->try_emplace({k0,k1}, fwd1D.at(j).at(i));
            }
            nodeId++;
         }
      }
   }

   void KaHIPSplitting::computeEdgeCut()
   {
      std::vector<int> xadj;
      std::vector<int> adjncy;
      std::vector<int> vwgt;
      std::vector<int> adjcwgt;
      std::vector<int> xnodes = {0};

      // create 1D nodes
      std::cerr << "INDEXES FOR 1D" << std::endl;
      std::map<std::pair<int,int>,std::vector<int>> map1Dgrid;
      this->mapNodes(this->mMap1Dnodes, nullptr, &map1Dgrid, Dimensions::Transform::TRA1D, xnodes.back());
      xnodes.push_back(xnodes.back() + this->mMap1Dnodes.size());

      // create 2D nodes
      std::cerr << "INDEXES FOR 2D" << std::endl;
      std::map<std::pair<int,int>,std::vector<int>> map2Dmodes;
      std::map<std::pair<int,int>,std::vector<int>> map2Dgrid;
      this->mapNodes(this->mMap2Dnodes, &map2Dmodes, &map2Dgrid, Dimensions::Transform::TRA2D, xnodes.back());
      xnodes.push_back(xnodes.back() + this->mMap2Dnodes.size());

      // create 3D nodes
      std::cerr << "INDEXES FOR 3D" << std::endl;
      std::map<std::pair<int,int>,std::vector<int>> map3Dmodes;
      this->mapNodes(this->mMap3Dnodes, &map3Dmodes, nullptr, Dimensions::Transform::TRA3D, xnodes.back());
      xnodes.push_back(xnodes.back() + this->mMap3Dnodes.size());

      bool weightVertex = false;
      int vweight = 1;
      int eweight = 1;

      // Connect nodes
      int sze = 0;
      xadj.push_back(sze);
      vwgt.push_back(1);
      for(auto&& [k, v]: mMap1Dnodes)
      {
         for(auto&& k2: map1Dgrid.at(k))
         {
            auto key = std::make_pair(k.second, k2);
            if(mMap2Dnodes.count(key) > 0)
            {
               adjncy.push_back(mMap2Dnodes.at(key));
               adjcwgt.push_back(eweight);
               sze++;
            }
            else
            {
               throw std::logic_error("Graph edges for 1D -> 2D are inconsistent");
            }
         }
         xadj.push_back(sze);

         // Vertex weight
         if(weightVertex)
         {
            vweight = 1;
         }
         else
         {
            vweight = 1;
         }
         vwgt.push_back(vweight);
      }
      for(auto&& [k, v]: mMap2Dnodes)
      {
         for(auto&& k2: map2Dmodes.at(k))
         {
            auto key = std::make_pair(k2, k.first);
            if(mMap1Dnodes.count(key) > 0)
            {
               adjncy.push_back(mMap1Dnodes.at(key));
               adjcwgt.push_back(eweight);
               sze++;
            }
            else
            {
               throw std::logic_error("Graph edges for 2D -> 1D are inconsistent");
            }
         }
         for(auto&& k2: map2Dgrid.at(k))
         {
            auto key = std::make_pair(k.second, k2);
            if(mMap3Dnodes.count(key) > 0)
            {
               adjncy.push_back(mMap3Dnodes.at(key));
               adjcwgt.push_back(eweight);
               sze++;
            }
            else
            {
               throw std::logic_error("Graph edges for 2D -> 3D are inconsistent");
            }
         }
         xadj.push_back(sze);

         // Vertex weight
         if(weightVertex)
         {
            vweight = 1;
         }
         else
         {
            vweight = 1;
         }
         vwgt.push_back(vweight);
      }
      for(auto&& [k, v]: mMap3Dnodes)
      {
         for(auto&& k2: map3Dmodes.at(k))
         {
            auto key = std::make_pair(k2, k.first);
            if(mMap2Dnodes.count(key) > 0)
            {
               adjncy.push_back(mMap2Dnodes.at(key));
               adjcwgt.push_back(eweight);
               sze++;
            }
            //else
            //{
            //   throw std::logic_error("Graph edges for 3D -> 2D are inconsistent");
            //}
         }
         xadj.push_back(sze);

         // Vertex weight
         if(weightVertex)
         {
            vweight = 1;
         }
         else
         {
            vweight = 1;
         }
         vwgt.push_back(vweight);
      }

      std::string filebase = "graph";
      details::writeMetis(filebase + ".metis", xadj, adjncy, vwgt, adjcwgt);

      int n            = xadj.size()-1;
      double imbalance = 0.03;
      mPartition.resize(n);
      int edge_cut     = 0;
      int nparts       = this->nCpu();
      kaffpa_balance_NE(&n, vwgt.data(), xadj.data(), adjcwgt.data(), adjncy.data(), &nparts, &imbalance, false, 0, STRONG, & edge_cut, mPartition.data());
      details::writePartition(filebase + "_partition.txt", mPartition);

      details::writeDot(filebase + ".dot", xnodes, xadj, adjncy, mPartition);
      details::writeDot(filebase + "_1D2D.dot", xnodes, xadj, adjncy, mPartition, 0);
      details::writeDot(filebase + "_2D3D.dot", xnodes, xadj, adjncy, mPartition, 1);
   }

   SharedTransformResolution  KaHIPSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status)
   {
      std::cerr << "SPLITTING DIMENSIONS " << static_cast<int>(transId) << " WITH KAHIP: rank = " << cpuId << std::endl;
      if(this->mPartition.size() == 0)
      {
         this->computeEdgeCut();
      }

      // Storage for the forward 1D indexes
      std::vector<std::vector<std::vector<int> > >  fwd1D;
      // Storage for the backward 1D indexes
      std::vector<std::vector<std::vector<int> > >  bwd1D;
      // Storage for the 2D indexes
      std::vector<std::vector<int> >  idx2D;
      // Storage for the 3D indexes
      std::vector<int>  idx3D;

      std::map<std::pair<int,int>,int> *pmapNodes;
      int i = 0;
      if(transId == Dimensions::Transform::TRA1D || transId == Dimensions::Transform::SPECTRAL)
      {
         std::cerr << "SPLIT TRA1D: rank = " << cpuId << std::endl;
         pmapNodes = &mMap1Dnodes;
         i = 0;
      }
      else if(transId == Dimensions::Transform::TRA2D)
      {
         std::cerr << "SPLIT TRA2D: rank = " << cpuId << std::endl;
         pmapNodes = &mMap2Dnodes;
         i = this->mMap1Dnodes.size();
      }
      else if(transId == Dimensions::Transform::TRA3D)
      {
         std::cerr << "SPLIT TRA3D: rank = " << cpuId << std::endl;
         pmapNodes = &mMap3Dnodes;
         i = this->mMap1Dnodes.size() + this->mMap2Dnodes.size();
      }
      std::map<int, std::vector<int>> filter3D;
      for(auto&& [k, v]: *pmapNodes)
      {
         if(this->mPartition.at(i) == cpuId)
         {
            if(filter3D.count(k.first) == 0)
            {
               filter3D.try_emplace(k.first, std::vector<int>());
            }
            filter3D.at(k.first).push_back(k.second);
         }
         i++;
      }
      for(auto&& [k, v]: filter3D)
      {
         idx3D.push_back(k);
         idx2D.push_back(v);
      }

      this->mspScheme->fillIndexes1D(transId, fwd1D, bwd1D, idx2D, idx3D);

      for(int j = 0; j < idx3D.size(); j++)
      {
         std::cerr << idx3D.at(j) << ": ";
         for(auto&& k: idx2D.at(j))
         {
            std::cerr << k << " ";
         }
      std::cerr << std::endl;
      }
      std::cerr << "%%%%%%%%%%%%%%%%%%%%%%%% DONE $$$$$$$$$$$$$$$$$$$" << std::endl;

      // Create TransformResolution object
      auto spTraRes = std::make_shared<TransformResolution>(fwd1D, bwd1D, idx2D, idx3D);
      std::cerr << "================= Resolution ===============" << std::endl;
      return spTraRes;
   }

   void KaHIPSplitting::selectGrouper(const Splitting::Groupers::Id selected)
   {
      // Only split in first transpose
      if(this->factors()(1) == 1 && selected != Splitting::Groupers::SINGLE1D)
      {
         this->mGrouper = Splitting::Groupers::EQUATION;
      }
      // Only split in second transpose
      else if(this->factors()(0) == 1 && selected != Splitting::Groupers::SINGLE2D)
      {
         this->mGrouper = Splitting::Groupers::EQUATION;
      }
      else
      {
         this->mGrouper = selected;
      }
   }

   Array KaHIPSplitting::computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp)
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

      std::cerr << comm.transpose() << std::endl;
      std::cerr << balance.transpose() << std::endl;
      std::cerr << details.transpose() << std::endl;

      return details;
   }

namespace details
{
void writeDot(std::string filename, const std::vector<int>& xnodes, const std::vector<int>& xadj, const std::vector<int>& adjncy, const std::vector<int>& part, const int stage)
{
   int n = xadj.size()-1;

   std::ofstream dot(filename);

   dot << "strict graph {" << std::endl;
   dot << "node [colorscheme=set19]" << std::endl;
   for(int s = 0; s < xnodes.size() - 1; s++)
   {
      for(int node = 0; node < xnodes[s+1] - xnodes[s]; node++)
      {
         dot << "\"s" << s << "_" << node << "\"" << " [style = filled, color=" << part[node] + 1 << "]" << std::endl;
      }
   }
   for(int node = 0; node < n; node++)
   {
      for(int i = xadj[node]; i < xadj[node+1]; i++)
      {
         int nLeft = node;
         int nRight = adjncy[i];
         int sLeft = -42;
         int sRight = -42;
         for(int k = 1; k < xnodes.size(); k++)
         {
            if(sLeft < 0 && nLeft - xnodes[k] < 0)
            {
               sLeft = k-1;
               nLeft -= xnodes[k-1];
            }
            if(sRight < 0 && nRight - xnodes[k] < 0)
            {
               sRight = k-1;
               nRight -= xnodes[k-1];
            }
         }

         if(stage < 0 || ((sLeft == stage && sRight == stage + 1) || (sLeft == stage + 1 && sRight == stage)))
         {
            dot << "\"s" << sLeft << "_" << nLeft << "\" -- \"s" << sRight << "_" << nRight << "\"" << std::endl;
         }
      }
   }

  // Close the file
  dot << "}" << std::endl;
  dot.close();
}

void writePartition(const std::string filename, const std::vector<int>& part)
{
   std::ofstream partition(filename);

   for(auto&& c: part)
   {
      partition << c << std::endl;
   }

   partition.close();
}

void writeMetis(const std::string filename, const std::vector<int>& xadj, const std::vector<int>& adjncy, const std::vector<int>& vwgt, const std::vector<int>& adjcwgt)
{
   std::ofstream metis(filename);

   metis << xadj.size() - 1 << " " << adjncy.size()/2 << " " << 11 << std::endl;

   for(int node = 0; node < xadj.size()-1; node++)
   {
      metis << vwgt[node];
      for(int i = xadj[node]; i < xadj[node+1]; i++)
      {
         metis << " " << adjncy[i] + 1 << " " << adjcwgt[i];
      }
      metis << std::endl;
   }

   metis.close();
}
}

}
}
