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
#include "QuICC/TestSuite/Framework/Communicators/TestHelper.hpp"
#include "QuICC/TestSuite/Framework/Communicators/TestArgs.hpp"
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
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Communicators {

   ArrayI processCmdLine()
   {
      // Set default arguments if required
      if(args().useDefault)
      {
         args().dim1D = 15;
         args().dim2D = 31;
         args().dim3D = 31;

         if(args().algorithm == "")
         {
            args().algorithm = "tubular";
         }

         if(args().grouper == "")
         {
            args().grouper = "transform";
         }

         args().params.clear();
      }
      else
      {
         if(args().algorithm == "")
         {
            args().algorithm = "serial";
         }

         if(args().grouper == "")
         {
            args().grouper = "equation";
         }
      }

      if(args().dim1D == 0 || args().dim2D == 0 || args().dim3D == 0)
      {
         throw std::logic_error("Dimensions are not set properly");
      }

      // Set simulation truncation
      QuICC::ArrayI dim(3);
      dim << args().dim1D, args().dim2D, args().dim3D;

      INFO( "Input parameters" );
      INFO( "dim1D: " << dim(0) );
      INFO( "dim2D: " << dim(1) );
      INFO( "dim3D: " << dim(2) );

      return dim;
   }

   void initCommunicator(Test& test, const Parallel::SplittingDescription& descr)
   {
      auto spRes = test.spRes;
      auto&& comm = test.comm;
      auto&& fwdTree = test.fwdTree;
      auto&& bwdTree = test.bwdTree;
      auto spBwdGrouper = test.spBwdGrouper;
      auto spFwdGrouper = test.spFwdGrouper;

      std::vector<Dimensions::Transform::Id> tIds = {Dimensions::Transform::TRA1D, Dimensions::Transform::TRA2D, Dimensions::Transform::TRA3D};
      for(auto t: tIds)
      {
         comm.init(t, spRes->spFwdSetup(t), spRes->spBwdSetup(t));
      }
      comm.init(Dimensions::Transform::SPECTRAL, spRes->spSpectralSetup(), spRes->spSpectralSetup());

      std::vector<Transform::TransformPath> t;
      std::map<size_t, std::vector<Transform::TransformPath> > mt;

      auto spSteps = Transform::createTransformSteps(spRes->sim().spSpatialScheme());
      // Create forward scalar transform tree
      std::vector<Transform::ITransformSteps::PathId> comps = {{FieldComponents::Spectral::SCALAR,Transform::Path::Scalar::id()}};
      t = spSteps->forwardScalar(comps);
      mt.insert(std::make_pair(PhysicalNames::Temperature::id(), t));

      Transform::TransformTreeTools::generateTrees(fwdTree, mt, TransformDirection::FORWARD, "scalar");
      t.clear();
      mt.clear();

      // Create backward scalar transform tree
      std::map<FieldComponents::Physical::Id,bool> req = {{FieldComponents::Physical::SCALAR,true}};
      t = spSteps->backwardScalar(req);
      mt.insert(std::make_pair(PhysicalNames::Temperature::id(), t));

      Transform::TransformTreeTools::generateTrees(bwdTree, mt, TransformDirection::BACKWARD, "scalar");
      t.clear();
      mt.clear();

      // Initialize converters
      std::vector<ArrayI> packs;
      Parallel::setGrouper(descr, spFwdGrouper, spBwdGrouper);
      // Get the buffer pack sizes for first dimension
      packs.push_back(spFwdGrouper->packs1D(fwdTree));
      packs.push_back(spBwdGrouper->packs1D(bwdTree));

      if(spRes->sim().ss().dimension() == 3)
      {
         // Get the buffer pack sizes for second dimension
         packs.push_back(spFwdGrouper->packs2D(fwdTree));
         packs.push_back(spBwdGrouper->packs2D(bwdTree));
      }
      comm.initConverter(spRes, packs, spFwdGrouper->split);
   }

   void setupSpectralCommunication(Test& test)
   {
      auto&& comm = test.comm;

      comm.template converter<QuICC::Dimensions::Transform::TRA1D>().setupCommunication(test.spBwdGrouper->packs1D(test.bwdTree.at(0)), QuICC::TransformDirection::BACKWARD);
      comm.template converter<QuICC::Dimensions::Transform::TRA1D>().prepareBackwardReceive();
   }

   void setup1D2D3DCommunication(Test& test)
   {
      auto&& comm = test.comm;
      auto&& bwdTree = test.bwdTree;
      auto spBwdGrouper = test.spBwdGrouper;

      CHECK( bwdTree.size() == 1 );
      comm.template converter<QuICC::Dimensions::Transform::TRA2D>().setupCommunication(spBwdGrouper->packs1D(bwdTree.at(0)), QuICC::TransformDirection::BACKWARD);
      comm.template converter<QuICC::Dimensions::Transform::TRA2D>().prepareBackwardReceive();
      comm.template converter<QuICC::Dimensions::Transform::TRA3D>().setupCommunication(spBwdGrouper->packs2D(bwdTree.at(0)), QuICC::TransformDirection::BACKWARD);
      comm.template converter<QuICC::Dimensions::Transform::TRA3D>().prepareBackwardReceive();
      }

   void transposeSpectral_1D(Test& test)
   {
      auto&& res = *test.spRes;
      auto&& comm = test.comm;

      INFO( "Checking spectral -> 1D transpose" );
      const auto inId = Dimensions::Transform::SPECTRAL;
      auto pOutData = res.sim().ss().bwdPtr(inId);
      comm.storage(inId).provideFwd(pOutData);
      MHDFloat scale = std::pow(10,1+std::floor(std::log10(res.sim().dimensions(Dimensions::Space::SPECTRAL).maxCoeff())));
      INFO( "Scale: " << scale );
      // Initialize pOutData
      std::visit(
            [&](auto&& p)
            {
               const auto& tRes = *res.cpu()->dim(inId);
               int *si, *sj, *sk;
               int i_, j_, k_;
               // Extract radial dimension for ijk
               if(res.sim().ss().has(SpatialScheme::Feature::SpectralOrdering123))
               {
                  si = &i_;
                  sj = &j_;
                  sk = &k_;
               }
               // Extract radial dimension for ikj
               else
               {
                  si = &i_;
                  sj = &k_;
                  sk = &j_;
               }
               // Initialize to bad value
               p->rData().setConstant(badValueOut);
               // Initialize with id representing mode, c depends on resolution: id = (i + 1) + c*j + c^2*k
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
                     for(int i = 0; i < tRes.dim<Dimensions::Data::DATB1D>(k); i++)
                     {
                        i_ = tRes.idx<Dimensions::Data::DATB1D>(i,k);
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
      comm.converter<Dimensions::Transform::TRA1D>().initiateForwardSend();
      comm.receiveBackward(outId, pInData);

      // Check received data
      std::visit(
            [&](auto&& p)
            {
               int *si, *sj, *sk;
               int i_, j_, k_;
               if(res.sim().ss().has(SpatialScheme::Feature::TransformSpectralOrdering123))
               {
                  si = &i_;
                  sj = &j_;
                  sk = &k_;
               }
               else
               {
                  si = &i_;
                  sj = &k_;
                  sk = &j_;
               }
               const auto& tRes = *res.cpu()->dim(outId);
               const auto& cnt = res.counter();
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
                     CHECK( tRes.dim<Dimensions::Data::DATB1D>(k) >= cnt.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, *sj) );
                     for(int i = 0; i < cnt.dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL, *sj); i++)
                     {
                        i_ = tRes.idx<Dimensions::Data::DATB1D>(i,k);
                        MHDFloat modeId = (*si) + 1 + (*sj)*scale + (*sk)*scale*scale;
                        INFO( "i,j,k: " << i << "," << j << "," << k );
                        INFO( "si,sj,sk: " << *si << "," << *sj << "," << *sk );
                        CHECK( std::abs(p->point(i,j,k)) == modeId );
                     }
                  }
               }
            },
         pInData);
   }

   void transpose1D_2D(Test& test)
   {
      auto&& res = *test.spRes;
      auto&& comm = test.comm;

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
               const auto& tRes = *res.cpu()->dim(inId);
               int *si, *sj, *sk;
               int i_, j_, k_;
               if(res.sim().ss().has(SpatialScheme::Feature::TransformSpectralOrdering123))
               {
                  si = &i_;
                  sj = &j_;
                  sk = &k_;
               }
               else
               {
                  si = &i_;
                  sj = &k_;
                  sk = &j_;
               }
               // Initialize to bad value
               p->rData().setConstant(badValueOut);
               // Initialize with id representing mode, c depends on resolution: id = (i + 1) + c*j + c^2*k
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
                     for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(k); i++)
                     {
                        i_ = tRes.idx<Dimensions::Data::DATF1D>(i,k);
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
               const auto& tRes = *res.cpu()->dim(outId);
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
                     for(int i = 0; i < tRes.dim<Dimensions::Data::DATB1D>(k); i++)
                     {
                        i_ = tRes.idx<Dimensions::Data::DATB1D>(i,k);
                        MHDFloat modeId = (*si) + 1 + (*sj)*scale + (*sk)*scale*scale;
                        CHECK( std::abs(p->point(i,j,k)) == modeId );
                     }
                  }
               }
            },
         pInData);
   }

   void transpose2D_3D(Test& test)
   {
      auto&& res = *test.spRes;
      auto&& comm = test.comm;

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
               const auto& tRes = *res.cpu()->dim(inId);
               int *si, *sj, *sk;
               int i_, j_, k_;
               si = &j_;
               sj = &i_;
               sk = &k_;
               // Initialize to bad value
               p->rData().setConstant(badValueOut);
               // Initialize with id representing mode, c depends on resolution: id = (i + 1) + c*j + c^2*k
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);
                     for(int i = 0; i < tRes.dim<Dimensions::Data::DATF1D>(k); i++)
                     {
                        i_ = tRes.idx<Dimensions::Data::DATF1D>(i,k);
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
               const auto& tRes = *res.cpu()->dim(outId);
               for(int k = 0; k < tRes.dim<Dimensions::Data::DAT3D>(); k++)
               {
                  k_ = tRes.idx<Dimensions::Data::DAT3D>(k);
                  for(int j = 0; j < tRes.dim<Dimensions::Data::DAT2D>(k); j++)
                  {
                     j_ = tRes.idx<Dimensions::Data::DAT2D>(j,k);

                     // Check transfered modes
                     for(int i = 0; i < res.sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL); i++)
                     {
                        i_ = tRes.idx<Dimensions::Data::DATB1D>(i,k);
                        MHDFloat modeId = (*si) + 1 + (*sj)*scale + (*sk)*scale*scale;
                        CHECK( std::abs(p->point(i,j,k)) == modeId );
                     }

                     // Check dealiased modes
                     for(int i = res.sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL); i < tRes.dim<Dimensions::Data::DATB1D>(k); i++)
                     {
                        i_ = tRes.idx<Dimensions::Data::DATB1D>(i,k);
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
