/**
 * @file @Scheme@CommunicatorTest.cpp
 * @brief Tests for the communiator
 */

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>
#include <fstream>

// Project includes
//
#include "QuICC/Communicators/Communicator.hpp"
#include "QuICC/TestSuite/Framework/Communicators/CommunicatorStages.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/Feature.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/TransformConfigurators/TransformStepsFactory.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/TransformConfigurators/TransformTreeTools.hpp"
#include "QuICC/TypeSelectors/ParallelSelector.hpp"
#include "QuICC/Transform/Path/Scalar.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Communicators {

   void initCommunicator(CommCoordinator& coord, SharedResolution spRes, const Parallel::SplittingDescription& descr)
   {
      auto&& comm = coord.comm;

      std::vector<Dimensions::Transform::Id> tIds = {Dimensions::Transform::TRA1D, Dimensions::Transform::TRA2D, Dimensions::Transform::TRA3D};
      for(auto t: tIds)
      {
         comm.init(t, spRes->spFwdSetup(t), spRes->spBwdSetup(t));
      }

      std::vector<Transform::TransformPath> t;
      std::map<size_t, std::vector<Transform::TransformPath> > mt;

      auto spSteps = Transform::createTransformSteps(spRes->sim().spSpatialScheme());
      // Create forward scalar transform tree
      std::vector<Transform::ITransformSteps::PathId> comps = {{FieldComponents::Spectral::SCALAR,Transform::Path::Scalar::id()}};
      t = spSteps->forwardScalar(comps);
      mt.insert(std::make_pair(PhysicalNames::Temperature::id(), t));

      Transform::TransformTreeTools::generateTrees(coord.fwdTree, mt, TransformDirection::FORWARD, "scalar");
      t.clear();
      mt.clear();

      // Create backward scalar transform tree
      std::map<FieldComponents::Physical::Id,bool> req = {{FieldComponents::Physical::SCALAR,true}};
      t = spSteps->backwardScalar(req);
      mt.insert(std::make_pair(PhysicalNames::Temperature::id(), t));

      Transform::TransformTreeTools::generateTrees(coord.bwdTree, mt, TransformDirection::BACKWARD, "scalar");
      t.clear();
      mt.clear();

      // Initialize converters
      std::vector<ArrayI> packs;
      Parallel::setGrouper(descr, coord.spFwdGrouper, coord.spBwdGrouper);
      // Get the buffer pack sizes for first dimension
      packs.push_back(coord.spFwdGrouper->packs1D(coord.fwdTree));
      packs.push_back(coord.spBwdGrouper->packs1D(coord.bwdTree));

      if(spRes->sim().ss().dimension() == 3)
      {
         // Get the buffer pack sizes for second dimension
         packs.push_back(coord.spFwdGrouper->packs2D(coord.fwdTree));
         packs.push_back(coord.spBwdGrouper->packs2D(coord.bwdTree));
      }
      comm.initConverter(spRes, packs, coord.spFwdGrouper->split);
   }

   void transposeSpectral_1D(const Resolution& res, Parallel::Communicator& comm)
   {
      INFO( "Checking spectral -> 1D transpose" );
      const auto inId = Dimensions::Transform::TRA1D;
      auto pOutData = res.sim().ss().bwdPtr(inId);
      comm.storage(inId).provideBwd(pOutData);
      MHDFloat scale = std::pow(10,1+std::floor(std::log10(res.sim().dimensions(Dimensions::Space::SPECTRAL).maxCoeff())));
      INFO( "Scale: " << scale );
      // Initialize pOutData
      std::map<int,int> nN;
      std::visit(
            [&](auto&& p)
            {
            auto spTRes = res.cpu()->dim(inId);
            int *si, *sj, *sk;
            int i_, j_, k_;
            // Extract radial dimension for ijk
            if(res.sim().ss().has(SpatialScheme::Feature::SpectralOrdering123))
            {
               si = &i_;
               sj = &j_;
               sk = &k_;
               // This will fail with triangular truncation ??? size depends on l not m
               for(int k = 0; k < spTRes->dim<Dimensions::Data::DAT3D>(); k++)
               {
                  for(int j = 0; j < spTRes->dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = spTRes->idx<Dimensions::Data::DAT2D>(j,k);
                     nN.emplace(*sj, spTRes->dim<Dimensions::Data::DATB1D>(k));
                  }
               }
            // Extract radial dimension for ikj
            } else
            {
               si = &i_;
               sj = &k_;
               sk = &j_;
               for(int k = 0; k < spTRes->dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = spTRes->idx<Dimensions::Data::DAT3D>(k);
                  nN.emplace(*sj, spTRes->dim<Dimensions::Data::DATB1D>(k));
               }
            }
            // Initialize to bad value
            p->rData().setConstant(badValueOut);
            // Initialize with id representing mode, c depends on resolution: id = (i + 1) + c*j + c^2*k
            for(int k = 0; k < spTRes->dim<Dimensions::Data::DAT3D>(); k++)
            {
               k_ = spTRes->idx<Dimensions::Data::DAT3D>(k);
               for(int j = 0; j < spTRes->dim<Dimensions::Data::DAT2D>(k); j++)
               {
                  j_ = spTRes->idx<Dimensions::Data::DAT2D>(j,k);
                  for(int i = 0; i < spTRes->dim<Dimensions::Data::DATB1D>(k); i++)
                  {
                     i_ = spTRes->idx<Dimensions::Data::DATB1D>(i,k);
                     MHDFloat modeId = (*si) + 1 + (*sj)*scale + (*sk)*scale*scale;
                     // Safety check to make sure modeId is never 0
                     CHECK( modeId != 0 );
                     p->rPoint(i,j,k) = modeId;
                  }
               }
            }
            },
         pOutData);

      // Prepare input data
      const auto outId = Dimensions::Transform::TRA1D;
      auto pInData = res.sim().ss().bwdPtr(outId);

      // Transfer and receive
      comm.storage(inId).holdBwd(pOutData);
      comm.receiveBackward(outId, pInData);

      // Check received data
      std::visit(
            [&](auto&& p)
            {
               int *si, *sj, *sk;
               int i_, j_, k_;
               if(res.sim().ss().has(SpatialScheme::Feature::SpectralOrdering123))
               {
                  si = &i_;
                  sj = &j_;
                  sk = &k_;
               } else
               {
                  si = &i_;
                  sj = &k_;
                  sk = &j_;
               }
               auto spTRes = res.cpu()->dim(outId);
               for(int k = 0; k < spTRes->dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = spTRes->idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < spTRes->dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = spTRes->idx<Dimensions::Data::DAT2D>(j,k);
                     CHECK( spTRes->dim<Dimensions::Data::DATB1D>(k) >= nN.at(*sj) );
                     for(int i = 0; i < nN.at(*sj); i++)
                     {
                        i_ = spTRes->idx<Dimensions::Data::DATB1D>(i,k);
                        MHDFloat modeId = (*si) + 1 + (*sj)*scale + (*sk)*scale*scale;
                        CHECK( std::abs(p->point(i,j,k)) == modeId );
                     }
                  }
               }
            },
         pInData);
   }

   void transpose1D_2D(const Resolution& res, Parallel::Communicator& comm)
   {
      INFO( "Checking 1D -> 2D transpose" );
      const auto inId = Dimensions::Transform::TRA1D;
      auto pOutData = res.sim().ss().fwdPtr(inId);
      comm.storage(inId).provideFwd(pOutData);
      MHDFloat scale = std::pow(10,1+std::floor(std::log10(res.sim().dimensions(Dimensions::Space::SPECTRAL).maxCoeff())));
      INFO( "Scale: " << scale );
      // Initialize pOutData
      std::visit(
            [&](auto&& p)
            {
               auto spTRes = res.cpu()->dim(inId);
               int *si, *sj, *sk;
               int i_, j_, k_;
               if(res.sim().ss().has(SpatialScheme::Feature::SpectralOrdering123))
               {
                  si = &i_;
                  sj = &j_;
                  sk = &k_;
               } else
               {
                  si = &i_;
                  sj = &k_;
                  sk = &j_;
               }
               // Initialize to bad value
               p->rData().setConstant(badValueOut);
               // Initialize with id representing mode, c depends on resolution: id = (i + 1) + c*j + c^2*k
               for(int k = 0; k < spTRes->dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = spTRes->idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < spTRes->dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = spTRes->idx<Dimensions::Data::DAT2D>(j,k);
                     for(int i = 0; i < spTRes->dim<Dimensions::Data::DATF1D>(k); i++)
                     {
                        i_ = spTRes->idx<Dimensions::Data::DATF1D>(i,k);
                        MHDFloat modeId = (*si) + 1 + (*sj)*scale + (*sk)*scale*scale;
                        // Safety check to make sure modeId is never 0
                        CHECK( modeId != 0 );
                        p->rPoint(i,j,k) = modeId;
                     }
                  }
               }
            },
         pOutData);


      // Prepare input data
      const auto outId = Dimensions::Transform::TRA2D;
      auto pInData = res.sim().ss().bwdPtr(outId);
      std::size_t addr;
      comm.storage(outId).provideBwd(pInData);
      std::visit(
            [&](auto&& p)
            {
               // store buffer address
               addr = reinterpret_cast<std::size_t>(p->data().data());
               // Initialize buffer with bad values
               p->rData().setConstant(badValueIn);
            },
         pInData);
      comm.storage(outId).freeBwd(pInData);

      // Transfer and receive
      comm.transferForward(inId, pOutData);
      comm.converter<Dimensions::Transform::TRA2D>().initiateForwardSend();
      comm.receiveBackward(outId, pInData);

      // Check received data
      std::visit(
            [&](auto&& p)
            {
               // Check new data is stored in the same buffer
               CHECK( addr == reinterpret_cast<std::size_t>(p->data().data()) );

               int *si, *sj, *sk;
               int i_, j_, k_;
               si = &j_;
               sj = &i_;
               sk = &k_;
               auto spTRes = res.cpu()->dim(outId);
               for(int k = 0; k < spTRes->dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = spTRes->idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < spTRes->dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = spTRes->idx<Dimensions::Data::DAT2D>(j,k);
                     for(int i = 0; i < spTRes->dim<Dimensions::Data::DATB1D>(k); i++)
                     {
                        i_ = spTRes->idx<Dimensions::Data::DATB1D>(i,k);
                        MHDFloat modeId = (*si) + 1 + (*sj)*scale + (*sk)*scale*scale;
                        CHECK( std::abs(p->point(i,j,k)) == modeId );
                     }
                  }
               }
            },
         pInData);
   }

   void transpose2D_3D(const Resolution& res, Parallel::Communicator& comm)
   {
      INFO( "Checking 2D -> 3D transpose" );
      const auto inId = Dimensions::Transform::TRA2D;
      auto pOutData = res.sim().ss().fwdPtr(inId);
      comm.storage(inId).provideFwd(pOutData);
      MHDFloat scale = std::pow(10,1+std::floor(std::log10(2.0*res.sim().dimensions(Dimensions::Space::SPECTRAL).maxCoeff())));
      INFO( "Scale: " << scale );
      // Initialize pOutData
      std::visit(
            [&](auto&& p)
            {
               auto spTRes = res.cpu()->dim(inId);
               int *si, *sj, *sk;
               int i_, j_, k_;
               si = &j_;
               sj = &i_;
               sk = &k_;
               // Initialize to bad value
               p->rData().setConstant(badValueOut);
               // Initialize with id representing mode, c depends on resolution: id = (i + 1) + c*j + c^2*k
               for(int k = 0; k < spTRes->dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = spTRes->idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < spTRes->dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = spTRes->idx<Dimensions::Data::DAT2D>(j,k);
                     for(int i = 0; i < spTRes->dim<Dimensions::Data::DATF1D>(k); i++)
                     {
                        i_ = spTRes->idx<Dimensions::Data::DATF1D>(i,k);
                        MHDFloat modeId = (*si) + 1 + (*sj)*scale + (*sk)*scale*scale;
                        // Safety check to make sure modeId is never 0
                        CHECK( modeId != 0 );
                        p->rPoint(i,j,k) = modeId;
                     }
                  }
               }
            },
         pOutData);

      // Prepare input data
      const auto outId = Dimensions::Transform::TRA3D;
      auto pInData = res.sim().ss().bwdPtr(outId);
      std::size_t addr;
      comm.storage(outId).provideBwd(pInData);
      std::visit(
            [&](auto&& p)
            {
               // store buffer address
               addr = reinterpret_cast<std::size_t>(p->data().data());
               // Initialize buffer with bad values
               p->rData().setConstant(badValueIn);
            },
         pInData);
      comm.storage(outId).freeBwd(pInData);

      // Transfer and receive
      comm.transferForward(inId, pOutData);
      comm.converter<Dimensions::Transform::TRA3D>().initiateForwardSend();
      comm.receiveBackward(outId, pInData);

      // Check received data
      std::visit(
            [&](auto&& p)
            {
               // Check new data is stored in the same buffer
               CHECK( addr == reinterpret_cast<std::size_t>(p->data().data()) );

               int *si, *sj, *sk;
               int i_, j_, k_;
               si = &k_;
               sj = &j_;
               sk = &i_;
               auto spTRes = res.cpu()->dim(outId);
               for(int k = 0; k < spTRes->dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = spTRes->idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < spTRes->dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = spTRes->idx<Dimensions::Data::DAT2D>(j,k);

                     // Check transfered modes
                     for(int i = 0; i < res.sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL); i++)
                     {
                        i_ = spTRes->idx<Dimensions::Data::DATB1D>(i,k);
                        MHDFloat modeId = (*si) + 1 + (*sj)*scale + (*sk)*scale*scale;
                        CHECK( std::abs(p->point(i,j,k)) == modeId );
                     }

                     // Check dealiased modes
                     for(int i = res.sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL); i < spTRes->dim<Dimensions::Data::DATB1D>(k); i++)
                     {
                        i_ = spTRes->idx<Dimensions::Data::DATB1D>(i,k);
                        CHECK( std::abs(p->point(i,j,k)) == badValueIn );
                     }
                  }
               }
            },
         pInData);
   }

}
}
}
}
