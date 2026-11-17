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
#include "Profiler/Interface.hpp"
#include "QuICC/PhysicalNames/registerAll.hpp"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/Pseudospectral/Utils.hpp"
#include "View/View.hpp"
#include "View/ViewUtils.hpp"
#include "ViewOps/ViewMemoryUtils.hpp"


namespace QuICC {

namespace Pseudospectral {

void Coordinator::addGraph(const std::string& graphStr,
   const Graph::PhysicalParameters<MHDFloat>& physParams)
{
   mGraphStr = graphStr;
   mGraphPhysParams = physParams;
}

void Coordinator::processGraph(const std::string& graphStr,
   const Graph::PhysicalParameters<MHDFloat>& physParams)
{
   // get Dims from mspRes
   // std::uint32_t Nr = jwRes.dim<Dimensions::Data::DATF1D>();
   std::uint32_t Nr = mspRes->sim().dim(Dimensions::Simulation::SIM1D,
      Dimensions::Space::PHYSICAL);
   std::uint32_t N = mspRes->sim().dim(Dimensions::Simulation::SIM1D,
      Dimensions::Space::SPECTRAL);
   // std::uint32_t Ntheta = alRes.dim<Dimensions::Data::DATF1D>();
   std::uint32_t Ntheta = mspRes->sim().dim(Dimensions::Simulation::SIM2D,
      Dimensions::Space::PHYSICAL);
   std::uint32_t L = mspRes->sim().dim(Dimensions::Simulation::SIM2D,
      Dimensions::Space::SPECTRAL);
   // std::uint32_t Nphi = ftRes.dim<Dimensions::Data::DATF1D>();
   std::uint32_t Nphi = mspRes->sim().dim(Dimensions::Simulation::SIM3D,
      Dimensions::Space::PHYSICAL);
   std::uint32_t M = mspRes->sim().dim(Dimensions::Simulation::SIM3D,
      Dimensions::Space::SPECTRAL);

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

   // Map view for scalar variables
   for(auto& [fId, s]: mScalarVariables)
   {
      std::visit(
         [&](auto&& p)
         {
            mId2View[fId] = p->rDom(0).rPerturbation().rGlobalView();

            std::array<std::uint32_t, dim> dims{modsDims[0], modsDims[2],
            modsDims[1]};

            // Return dimensions
            // mlir has layer first
            dimRets.push_back({dims[2], dims[0], dims[1]});
            layRets.push_back(layOpt[2][0]);
         }, s);
   }

   // Map view for vector variables
   for(auto& [k, v]: mVectorVariables)
   {
      std::visit(
         [&](auto&& p)
         {
            std::vector<FieldComponents::Spectral::Id> comps = {FieldComponents::Spectral::TOR, FieldComponents::Spectral::POL};
            for(auto&& c: comps)
            {
               std::size_t hComp = hash_combine(k, c);
               mId2View[hComp] = p->rDom(0).rPerturbation().rComp(c).rGlobalView();

               std::array<std::uint32_t, dim> dims{modsDims[0], modsDims[2],
               modsDims[1]};

               // Return dimensions
               // mlir has layer first
               dimRets.push_back({dims[2], dims[0], dims[1]});
               layRets.push_back(layOpt[2][0]);
            }
         }, v);
   }

   if (mVectorEquations.at(/*it=*/0).size() == 2)
   {
      mIsMag = true;
   }

   // Physical space (aka FT space, Stage::PPP and Stage::MPP, QuICC Stage2)

#if 0
   for(auto& [hComp, s]: mScalarVariables)
   {
      std::visit(
         [&](auto&& p)
         {
               mId2View[hComp] = p->rDom(0).rPhys().rGlobalView();

               if(p->rDom(0).rPhys().rGlobalView().pointers()[1].size() != metaFT.ptr.size())
               {
                  throw std::logic_error("SETUP IN NOT CORRECT");
               }
               if(p->rDom(0).rPhys().rGlobalView().indices()[1].size() != metaFT.idx.size())
               {
                  throw std::logic_error("SETUP IN NOT CORRECT");
               }

               // Map view for scalar variables
               std::array<std::uint32_t, dim> dims{physDims[2], physDims[1],
               physDims[0]};

               // Return dimensions
               // mlir has layer first
               dimRets.push_back({dims[2], dims[0], dims[1]});
               layRets.push_back(layOpt[0][0]);
         }, s);
   }
#endif

   // Map view for vector variables
   for(auto& [k, v]: mVectorVariables)
   {
      std::visit(
         [&](auto&& p)
         {
            std::vector<FieldComponents::Physical::Id> comps = {FieldComponents::Physical::R, FieldComponents::Physical::THETA, FieldComponents::Physical::PHI};
            for(auto&& c: comps)
            {
               std::size_t hComp = hash_combine(k, c);
               mId2View[hComp] = p->rDom(0).rPhys().rComp(c).rGlobalView();

               // Map view for scalar variables
               std::array<std::uint32_t, dim> dims{physDims[2], physDims[1],
               physDims[0]};

               // Return dimensions
               // mlir has layer first
               dimRets.push_back({dims[2], dims[0], dims[1]});
               layRets.push_back(layOpt[0][0]);
            }
         }, v);
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
   opt.grouping.group = mGraphOptions.groupingSize;
   mJitter = std::make_unique<QuICC::Graph::Jit<3>>(graphStr, mMemRsr, physDims,
      modsDims, layOpt, Graph::Stage::MMM, Graph::Stage::MMM, meta, physParams,
      opt);
};

} // namespace Pseudospectral
} // namespace QuICC
