/**
 * @file MetisSplitting.cpp
 * @brief Source of the implementation of a load splitting algorithm using Metis
 */

// System includes
//
#include "metis.h"
#include <iostream>
#include <set>

// Project includes
//
#include "QuICC/LoadSplitter/Algorithms/MetisSplitting.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/LoadSplitter/Algorithms/SplittingTools.hpp"

namespace QuICC {

namespace Parallel {

   MetisSplitting::MetisSplitting(const int id, const int nCpu, const ArrayI& dim, Splitting::Algorithms::Id algorithm, const std::list<int>& factors)
      : SplittingAlgorithm(id, nCpu, dim, algorithm)
   {
      // Initialise the NCpu factors
      this->initFactors(1);
   }

   bool MetisSplitting::applicable() const
   {
      bool status = true;

      // Check that all three dimensions are splittable by factors

      return status;
   }
   
   void MetisSplitting::mapNodes(std::map<std::pair<int,int>,int>& nodes, std::map<std::pair<int,int>,std::vector<int>>* pmodes, std::map<std::pair<int,int>,std::vector<int>>* pgrid, const Dimensions::Transform::Id transId, const int start)
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

   void MetisSplitting::computeEdgeCut()
   {
      std::vector<idx_t> xadj;
      std::vector<idx_t> adjncy;
      std::vector<idx_t> vwgt;
      std::vector<idx_t> adjcwgt;
      std::vector<idx_t> vsize;
      std::vector<idx_t> xnodes = {0};

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

      int vweight = 1;
      int vcomm = 1;
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

         vwgt.push_back(vweight);
         vwgt.push_back(1);
         vwgt.push_back(0);
         vwgt.push_back(0);
         vsize.push_back(vcomm);
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
         vwgt.push_back(vweight);
         vwgt.push_back(0);
         vwgt.push_back(1);
         vwgt.push_back(0);
         vsize.push_back(vcomm);
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

         vwgt.push_back(vweight);
         vwgt.push_back(0);
         vwgt.push_back(0);
         vwgt.push_back(1);
         vsize.push_back(vcomm);
      }

      std::string filebase = "graph";
      details::writeMetis(filebase + ".metis", xadj, adjncy, vwgt, adjcwgt);

      idx_t nvtxs = xadj.size()-1;
      idx_t ncon = 4;
      std::vector<real_t> ubvec = {1.001, 1.001, 1.001, 1.001};
      
      mPartition.resize(nvtxs);
      idx_t edge_cut     = 0;
      idx_t nparts       = this->nCpu();
      std::vector<idx_t> partition(nvtxs);

      std::vector<idx_t> options(METIS_NOPTIONS);
      int status = METIS_SetDefaultOptions(options.data());
      //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
      options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

      //status = METIS_PartGraphRecursive(&nvtxs, &ncon, xadj.data(), adjncy.data(), vwgt.data(), vsize.data(), adjcwgt.data(), &nparts, NULL, ubvec.data(), NULL, &edge_cut, partition.data());
      status = METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(), vwgt.data(), vsize.data(), adjcwgt.data(), &nparts, NULL, ubvec.data(), NULL, &edge_cut, partition.data());

      mPartition = partition;
      details::writePartition(filebase + "_partition.txt", mPartition);
      details::writeDot(filebase + ".dot", xnodes, xadj, adjncy, mPartition);
      details::writeDot(filebase + "_1D2D.dot", xnodes, xadj, adjncy, mPartition, 0);
      details::writeDot(filebase + "_2D3D.dot", xnodes, xadj, adjncy, mPartition, 1);
   }

   SharedTransformResolution  MetisSplitting::splitDimension(const Dimensions::Transform::Id transId, const int cpuId, int& status)
   {
      std::cerr << "SPLITTING DIMENSIONS " << static_cast<int>(transId) << " WITH METIS: rank = " << cpuId << std::endl;
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

   void MetisSplitting::selectGrouper(const Splitting::Groupers::Id selected)
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

   Array MetisSplitting::computeScore(SharedResolution spResolution, const Splitting::Groupers::Id grp)
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

      std::cerr << "comm: " << comm.transpose() << std::endl;
      std::cerr << "balance: " << balance.transpose() << std::endl;
      std::cerr << "details: " << details.transpose() << std::endl;

      return details;
   }

}
}
