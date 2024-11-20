/**
 * @file Coordinator.cpp
 * @brief Source of the high level pseudospectral coordinator
 */

// System includes
//
#include <algorithm>
#include <stdexcept>
#include <type_traits>

// Project includes
//
// #include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/Pseudospectral/Utils.hpp"
#include "QuICC/PhysicalNames/registerAll.hpp"
#include "View/View.hpp"
#include "View/ViewUtils.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"
#include "Profiler/Interface.hpp"


namespace QuICC {

namespace Pseudospectral {

void Coordinator::addGraph(const std::string& graphStr, const Graph::PhysicalParameters<MHDFloat>& physParams)
   {
      // get Dims from mspRes
      // std::uint32_t Nr = jwRes.dim<Dimensions::Data::DATF1D>();
      std::uint32_t Nr = mspRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::PHYSICAL);
      std::uint32_t N = mspRes->sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // std::uint32_t Ntheta = alRes.dim<Dimensions::Data::DATF1D>();
      std::uint32_t Ntheta = mspRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::PHYSICAL);
      std::uint32_t L = mspRes->sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
      // std::uint32_t Nphi = ftRes.dim<Dimensions::Data::DATF1D>();
      std::uint32_t Nphi = mspRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::PHYSICAL);
      std::uint32_t M = mspRes->sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Memory resource, depends on backend
      #ifdef QUICC_HAS_CUDA_BACKEND
      mMemRsr = std::make_shared<QuICC::Memory::Cuda::Malloc>();
      #else
      mMemRsr = std::make_shared<QuICC::Memory::Cpu::NewDelete>();
      #endif

      // get meta from mspRes
      const auto& jwRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA1D);
      auto metaJW = details::getMeta(jwRes, L, mMemRsr);
      const auto& alRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA2D);
      auto metaAL = details::getMeta(alRes, M, mMemRsr);
      const auto& ftRes = *mspRes->cpu()->dim(Dimensions::Transform::TRA3D);
      auto metaFT = details::getMeta(ftRes, Nr, mMemRsr);

      constexpr std::uint32_t dim = 3;
      /// RThetaPhi - v012
      std::array<std::uint32_t, dim> physDims{Nr, Ntheta, Nphi};
      /// @brief Spectral dimensions
      /// NLM - v012
      std::array<std::uint32_t, dim> modsDims{N, L, M};

      // Layouts, depends on backend
      std::array<std::array<std::string, 2>, 3> layOpt;
      #ifdef QUICC_HAS_CUDA_BACKEND
      layOpt[0] = {"DCCSC3D", "DCCSC3D"};
      layOpt[1] = {"DCCSC3DJIK", "S1CLCSC3DJIK"};
      layOpt[2] = {"DCCSC3DJIK", "DCCSC3DJIK"};
      #else
      layOpt[0] = {"DCCSC3D", "DCCSC3D"};
      layOpt[1] = {"DCCSC3D", "S1CLCSC3D"};
      layOpt[2] = {"DCCSC3D", "DCCSC3D"};
      #endif

      // Store meta stages to pass to Jitter
      std::vector<QuICC::View::ViewBase<std::uint32_t>> meta;
      meta.push_back({metaFT.ptr.data(), metaFT.ptr.size()});
      meta.push_back({metaFT.idx.data(), metaFT.idx.size()});
      meta.push_back({metaAL.ptr.data(), metaAL.ptr.size()});
      meta.push_back({metaAL.idx.data(), metaAL.idx.size()});
      meta.push_back({metaJW.ptr.data(), metaJW.ptr.size()});
      meta.push_back({metaJW.idx.data(), metaJW.idx.size()});

      // Wrapper pass options
      std::vector<std::vector<std::int64_t>> dimRets;
      std::vector<std::string> layRets;

      // Modal space (aka JW space, Stage::PMM and Stage::MMM, QuICC Stage0)
      std::array<View::ViewBase<std::uint32_t>, dim> pointersMods;
      pointersMods[1] = View::ViewBase<std::uint32_t>(metaJW.ptr.data(), metaJW.ptr.size());
      std::array<View::ViewBase<std::uint32_t>, dim> indicesMods;
      indicesMods[1] = View::ViewBase<std::uint32_t>(metaJW.idx.data(), metaJW.idx.size());

      // View for outputs/inputs
      std::size_t hVelTor = hash_combine(PhysicalNames::Velocity::id(),        FieldComponents::Spectral::TOR);
      std::size_t hVelPol = hash_combine(PhysicalNames::Velocity::id(),        FieldComponents::Spectral::POL);
      std::vector<size_t> fields = {PhysicalNames::Temperature::id(), hVelTor, hVelPol};

      if (mVectorEquations.at(/*it=*/0).size() == 2)
      {
         mIsMag = true;
      }

      if (mIsMag)
      {
         std::size_t hMagTor = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Spectral::TOR);
         fields.push_back(hMagTor);
         std::size_t hMagPol = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Spectral::POL);
         fields.push_back(hMagPol);
      }

      // Add Views and Storage for each component
      auto scalarVarPtr = mspRes->sim().ss().bwdPtr(Dimensions::Transform::SPECTRAL);

      std::visit(
         [&](auto&& p)
         {
            // Get field scalar type
            using fld_t = typename std::remove_reference_t<decltype(p->data())>::Scalar;

            for (size_t f = 0; f < fields.size(); ++f)
            {
               // Field Id
               auto fId = fields[f];

               // mem block
               Memory::MemBlock<fld_t> block(modsDims[0]*metaJW.idx.size(), mMemRsr.get());

               // view
               // for now this works only for JW space
               #ifdef QUICC_HAS_CUDA_BACKEND
               using jwLay_t = View::DCCSC3DJIK;
               #else
               using jwLay_t = View::DCCSC3D;
               #endif
               std::array<std::uint32_t, dim> dims {modsDims[0], modsDims[2], modsDims[1]};
               View::View<fld_t, jwLay_t> view(block.data(), block.size(), dims.data(), pointersMods.data(), indicesMods.data());

               // Store block
               mBlocksData.push_back(std::move(block));

               // Store view
               mId2View[fId] = view;

               // Return dimensions
               // mlir has layer first
               dimRets.push_back({dims[2], dims[0], dims[1]});
               layRets.push_back(layOpt[2][0]);
            }
         }, scalarVarPtr);

      // Physical space (aka FT space, Stage::PPP and Stage::MPP, QuICC Stage2)
      std::array<View::ViewBase<std::uint32_t>, dim> pointersPhys;
      pointersPhys[1] = View::ViewBase<std::uint32_t>(metaFT.ptr.data(), metaFT.ptr.size());
      std::array<View::ViewBase<std::uint32_t>, dim> indicesPhys;
      indicesPhys[1] = View::ViewBase<std::uint32_t>(metaFT.idx.data(), metaFT.idx.size());

      // Add Views for physical space
      /// \todo move cfl computation into graph and
      /// don't store physical space

      std::size_t hVelR = hash_combine(PhysicalNames::Velocity::id(), FieldComponents::Physical::R);
      std::size_t hVelTheta = hash_combine(PhysicalNames::Velocity::id(), FieldComponents::Physical::THETA);
      std::size_t hVelPhi = hash_combine(PhysicalNames::Velocity::id(), FieldComponents::Physical::PHI);
      std::vector<size_t> physFields = {hVelR, hVelTheta, hVelPhi};

      if (mIsMag)
      {
         std::size_t hMagR = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Physical::R);
         std::size_t hMagTheta = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Physical::THETA);
         std::size_t hMagPhi = hash_combine(PhysicalNames::Magnetic::id(), FieldComponents::Physical::PHI);
         physFields.push_back(hMagR);
         physFields.push_back(hMagTheta);
         physFields.push_back(hMagPhi);
      }

      for (size_t f = 0; f < physFields.size(); ++f)
      {
         using fld_t = MHDFloat;

         // Field Id
         auto fId = physFields[f];

         // mem block
         Memory::MemBlock<fld_t> block(physDims[2]*indicesPhys[1].size(), mMemRsr.get());

         // view
         // for now this works only for FT space
         std::array<std::uint32_t, dim> dims {physDims[2], physDims[1], physDims[0]};
         View::View<fld_t, View::DCCSC3D> view(block.data(), block.size(), dims.data(), pointersPhys.data(), indicesPhys.data());

         // Store block
         mBlocksData.push_back(std::move(block));

         // Store view
         mId2View[fId] = view;

         // Return dimensions
         // mlir has layer first
         dimRets.push_back({dims[2], dims[0], dims[1]});
         layRets.push_back(layOpt[0][0]);
      }

      // Store meta blocks
      mBlocksMeta.push_back(std::move(metaFT.ptr));
      mBlocksMeta.push_back(std::move(metaFT.idx));
      mBlocksMeta.push_back(std::move(metaAL.ptr));
      mBlocksMeta.push_back(std::move(metaAL.idx));
      mBlocksMeta.push_back(std::move(metaJW.ptr));
      mBlocksMeta.push_back(std::move(metaJW.idx));

      // Jitter
      Graph::PipelineOptions opt;
      opt.wrap.dimRets = dimRets;
      opt.wrap.layRets = layRets;
      mJitter = std::make_unique<QuICC::Graph::Jit<3>>(graphStr, mMemRsr, physDims, modsDims, layOpt, Graph::Stage::MMM, Graph::Stage::MMM, meta, physParams, opt);
   };

} // Pseudospectral
} // QuICC
