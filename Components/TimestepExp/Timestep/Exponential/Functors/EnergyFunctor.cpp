/**
 * @file EnergyR2viewCpu_t.cpp
 * @brief Wrapper of the Worland EnergyR2 Reductor
 */

// System includes
//

// Project includes
//
#include "Timestep/Exponential/Functors/EnergyFunctor.hpp"
#include "ViewOps/Worland/OpsBuilder.hpp"
#include "ViewOps/Worland/OpsQuadrature.hpp"
#include "View/ViewSparse.hpp"
#include "ViewOps/Pointwise/Cpu/Pointwise.hpp"
#include "ViewOps/Pointwise/Functors.hpp"
#include "ViewOps/Reduction/Cpu/Reduction.hpp"
#include "ViewOps/Quadrature/Op.hpp"
#include "ViewOps/Quadrature/Impl.hpp"
#include "ViewOps/Worland/Tags.hpp"
#include "ViewOps/Worland/Builder.hpp"
#include "Profiler/Interface.hpp"


#include <cstdio>
#include <iostream>
#include <memory>
namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

EnergyFunctor::EnergyFunctor(std::shared_ptr<Resolution> spRes, const int gSize, std::shared_ptr<Memory::memory_resource> mem)
   : spRes(spRes), mcGridSize(gSize), _mem(mem)
{}

void EnergyFunctor::init(mods_t in) const
{
   if(this->mPrjR2)
   {
      return;
   }

   // Energy calculation requires a different quadrature
   Internal::Array igrid, iweights;
   using namespace QuICC::Transform::Worland;
   // Check EnergyR2_t grid
   {
      OpsQuadrature<EnergyR2_t> quad;
      quad.computeQuadrature(igrid, iweights, mcGridSize);
   }
   // Check Energy_t grid
   {
      Internal::Array tgrid, tweights;
      OpsQuadrature<Energy_t> quad;
      quad.computeQuadrature(tgrid, tweights, mcGridSize);
      if(tgrid.size() > igrid.size())
      {
         igrid = tgrid;
         iweights = tweights;
      }
   }
   // Check Energy_t grid
   {
      Internal::Array tgrid, tweights;
      OpsQuadrature<EnergyD1R1_t> quad;
      quad.computeQuadrature(tgrid, tweights, mcGridSize);
      if(tgrid.size() > igrid.size())
      {
         igrid = tgrid;
         iweights = tweights;
      }
   }

    auto spSetup = spRes->spTransformSetup(Dimensions::Transform::TRA1D);
    std::uint32_t M = igrid.size();
    std::uint32_t K = in.dims()[0];

    std::uint32_t P = in.dims()[2];

    std::vector<std::uint32_t> layers;
    for(std::uint32_t p = 0; p < P; ++p)
    {
       if(in.pointers()[1][p + 1] > in.pointers()[1][p])
       {
          layers.push_back(p);
       }
    }

    constexpr size_t rank = 3;
    // dim 0 - Nr - radial points
    // dim 1 - R  - radial modes
    // dim 2 - L  - harmonic degree
    std::array<std::uint32_t, rank> dimensionsPrj {M, K, P};

    // Make projector view operators
    using namespace QuICC::Transform::Quadrature;
    using backendPrj_t = Cpu::ImplOp<phys_t, mods_t, op_t>;
    using derivedPrj_t = Op<phys_t, mods_t, op_t, backendPrj_t>;

    // Create energy R2 operator
    {
       mPrjR2 = std::make_unique<derivedPrj_t>(this->_mem);
       derivedPrj_t& derivedPrjOp = dynamic_cast<derivedPrj_t&>(*mPrjR2);
       derivedPrjOp.allocOp(dimensionsPrj, layers);
       auto prjView = derivedPrjOp.getOp();

       using namespace QuICC::Polynomial::Worland;
       OpsBuilder<op_t, EnergyR2_t, bwd_t> tBuilderBwd;
       tBuilderBwd.compute(prjView, igrid, Internal::Array());

       
       // dim 0 - R  - radial modes
       // dim 1 - Nr - radial points
       // dim 2 - L  - harmonic degree
       std::array<std::uint32_t, rank> dimensionsInt {K, M, P};

       // Make view integrator operator
       using backendInt_t = Cpu::ImplOp<mods_t, phys_t, op_t>;
       using derivedInt_t = Op<mods_t, phys_t, op_t, backendInt_t>;

       mIntR2 = std::make_unique<derivedInt_t>(this->_mem);
       derivedInt_t& derivedOp = dynamic_cast<derivedInt_t&>(*mIntR2);
       derivedOp.allocOp(dimensionsInt, layers);
       auto intView = derivedOp.getOp();

       OpsBuilder<op_t, EnergyR2_t, fwd_t> tBuilderFwd;
       tBuilderFwd.compute(intView, igrid, iweights);
    }

    // Create energy OverR1 operator
    {
       mPrjOverR1 = std::make_unique<derivedPrj_t>(this->_mem);
       derivedPrj_t& derivedPrjOp = dynamic_cast<derivedPrj_t&>(*mPrjOverR1);
       derivedPrjOp.allocOp(dimensionsPrj, layers);
       auto prjView = derivedPrjOp.getOp();

       using namespace QuICC::Polynomial::Worland;
       OpsBuilder<op_t, Energy_t, bwd_t> tBuilderBwd;
       tBuilderBwd.compute(prjView, igrid, Internal::Array());

       
       // dim 0 - R  - radial modes
       // dim 1 - Nr - radial points
       // dim 2 - L  - harmonic degree
       std::array<std::uint32_t, rank> dimensionsInt {K, M, P};

       // Make view integrator operator
       using backendInt_t = Cpu::ImplOp<mods_t, phys_t, op_t>;
       using derivedInt_t = Op<mods_t, phys_t, op_t, backendInt_t>;

       mIntOverR1 = std::make_unique<derivedInt_t>(this->_mem);
       derivedInt_t& derivedOp = dynamic_cast<derivedInt_t&>(*mIntOverR1);
       derivedOp.allocOp(dimensionsInt, layers);
       auto intView = derivedOp.getOp();

       OpsBuilder<op_t, Energy_t, fwd_t> tBuilderFwd;
       tBuilderFwd.compute(intView, igrid, iweights);
    }

    // Create energy OverR1D1R1 operator
    {
       mPrjOverR1D1R1 = std::make_unique<derivedPrj_t>(this->_mem);
       derivedPrj_t& derivedPrjOp = dynamic_cast<derivedPrj_t&>(*mPrjOverR1D1R1);
       derivedPrjOp.allocOp(dimensionsPrj, layers);
       auto prjView = derivedPrjOp.getOp();

       using namespace QuICC::Polynomial::Worland;
       OpsBuilder<op_t, EnergyD1R1_t, bwd_t> tBuilderBwd;
       tBuilderBwd.compute(prjView, igrid, Internal::Array());

       
       // dim 0 - R  - radial modes
       // dim 1 - Nr - radial points
       // dim 2 - L  - harmonic degree
       std::array<std::uint32_t, rank> dimensionsInt {K, M, P};

       // Make view integrator operator
       using backendInt_t = Cpu::ImplOp<mods_t, phys_t, op_t>;
       using derivedInt_t = Op<mods_t, phys_t, op_t, backendInt_t>;

       mIntOverR1D1R1 = std::make_unique<derivedInt_t>(this->_mem);
       derivedInt_t& derivedOp = dynamic_cast<derivedInt_t&>(*mIntOverR1D1R1);
       derivedOp.allocOp(dimensionsInt, layers);
       auto intView = derivedOp.getOp();

       OpsBuilder<op_t, EnergyD1R1_t, fwd_t> tBuilderFwd;
       tBuilderFwd.compute(intView, igrid, iweights);
    }

    // setup pointwise dot product
    mDot = std::make_unique<QuICC::Pointwise::Cpu::Op<QuICC::Pointwise::ComponentDotFunctor<double>, modsAbs_t, mods_t, mods_t>>(QuICC::Pointwise::ComponentDotFunctor<double>(1.0));

    std::uint32_t N = in.dims()[1];

    // dim 0 - R  - radial modes
    // dim 1 - M  - harmonic order
    // dim 2 - L  - harmonic degree
    std::array<std::uint32_t, rank> modsDims = {K, N, P};
    mModsData = std::move(Memory::MemBlock<std::complex<double>>(in.size(), _mem.get()));
    mEModsData = std::move(Memory::MemBlock<std::complex<double>>(in.size(), _mem.get()));

    // dim 0 - Nr - radial points
    // dim 1 - M  - harmonic order
    // dim 2 - L  - harmonic degree
    std::array<std::uint32_t, rank> physDims = {M, N, P};
    mPhysData = std::move(Memory::MemBlock<std::complex<double>>(M*in.size()/K, _mem.get()));

    using namespace QuICC::View;

    // set views
    mModsView = mods_t(mModsData.data(), mModsData.size(), modsDims.data(), in.pointers(), in.indices(), in.lds());
    mEModsView = mods_t(mEModsData.data(), mEModsData.size(), modsDims.data(), in.pointers(), in.indices(), in.lds());
    mPhysView = phys_t(mPhysData.data(), mPhysData.size(), physDims.data(), in.pointers(), in.indices(), in.lds());
    
    mModsAbsData = std::move(Memory::MemBlock<double>(in.size(), _mem.get()));
    mModsAbsView = modsAbs_t(mModsAbsData.data(), mModsAbsData.size(), modsDims.data(), in.pointers(), in.indices(), in.lds());

    // setup reductor
    mRed = std::make_unique<QuICC::Reduction::Cpu::Op<power_t, modsAbs_t, 0>>();

    // dim 0 - M  - harmonic order
    // dim 1 - L  - harmonic degree
    std::array<std::uint32_t, rank> modsRedDims = {N, P};

    // set pointers
    ViewBase<std::uint32_t> dataPointers2D[rank-1];
    dataPointers2D[0] = in.pointers()[1];
    ViewBase<std::uint32_t> dataIndices2D[rank-1];
    dataIndices2D[0] = in.indices()[1];

    // set views
    mModsRedData = std::move(Memory::MemBlock<double>(in.size()/K, _mem.get()));
    mModsRedView = power_t(mModsRedData.data(), mModsRedData.size(), modsRedDims.data(), dataPointers2D, dataIndices2D);
    
}

void EnergyFunctor::convert(mods_t out, mods_t in, const SpectralFieldId& id) const
{
   if(id.second == FieldComponents::Spectral::SCALAR || id.second == FieldComponents::Spectral::TOR || id.second == FieldComponents::Spectral::T)
   {
      /// project to physical space
      mPrjR2->apply(mPhysView, in);
      /// project to modal space
      mIntR2->apply(out, mPhysView);
   }
   else if(id.second == FieldComponents::Spectral::Q)
   {
      // Energy from q component
      /// project to physical space
      mPrjOverR1->apply(mPhysView, in);
      /// project to modal space
      mIntOverR1->apply(out, mPhysView);
   }
   else if(id.second == FieldComponents::Spectral::S)
   {
      // Energy from s component
      /// project to physical space
      mPrjOverR1D1R1->apply(mPhysView, in);
      /// project to modal space
      mIntOverR1D1R1->apply(out, mPhysView);
   }
   else
   {
      throw std::logic_error("Unknown spectral component");
   }
}

double EnergyFunctor::operator()(mods_t u, const SpectralFieldId& id) const
{
   double energy = 0.0;

   this->convert(mModsView, u, id);

   /// pointwise op on coefficients
   mDot->apply(mModsAbsView, mModsView, mModsView);
   /// reduction op on coefficients
   mRed->apply(mModsRedView, mModsAbsView);

   auto sM = std::make_unique<ScaleM>();
   if(id.second == FieldComponents::Spectral::SCALAR)
   {
      auto sL = std::make_unique<NoScale>();
      energy += this->computeEnergy(mModsRedView, sL, sM);
   }
   else if(id.second == FieldComponents::Spectral::Q)
   {
      auto sLq = std::make_unique<ScaleLq>();
      energy += this->computeEnergy(mModsRedView, sLq, sM);
   }
   else if(id.second == FieldComponents::Spectral::S)
   {
      auto sLs = std::make_unique<ScaleLst>();
      energy += this->computeEnergy(mModsRedView, sLs, sM);
   }
   else if(id.second == FieldComponents::Spectral::T)
   {
      auto sL = std::make_unique<ScaleLst>();
      energy += this->computeEnergy(mModsRedView, sL, sM);
   }
   else
   {
      throw std::logic_error("Unknown spectral component");
   }

   return energy;
}

double EnergyFunctor::operator()(mods_t u, mods_t vE, const SpectralFieldId& id) const
{
   double energy = 0.0;

   this->convert(mModsView, u, id);

   /// pointwise op on coefficients
   mDot->apply(mModsAbsView, mModsView, vE);
   /// reduction op on coefficients
   mRed->apply(mModsRedView, mModsAbsView);

   auto sM = std::make_unique<ScaleM>();
   if(id.second == FieldComponents::Spectral::SCALAR)
   {
      auto sL = std::make_unique<NoScale>();
      energy += this->computeEnergy(mModsRedView, sL, sM);
   }
   else if(id.second == FieldComponents::Spectral::Q)
   {
      auto sLq = std::make_unique<ScaleLq>();
      energy += this->computeEnergy(mModsRedView, sLq, sM);
   }
   else if(id.second == FieldComponents::Spectral::S)
   {
      auto sLs = std::make_unique<ScaleLst>();
      energy += this->computeEnergy(mModsRedView, sLs, sM);
   }
   else if(id.second == FieldComponents::Spectral::T)
   {
      auto sL = std::make_unique<ScaleLst>();
      energy += this->computeEnergy(mModsRedView, sL, sM);
   }
   else
   {
      throw std::logic_error("Unknown spectral component");
   }

   return energy;
}

double ScaleM::operator()(const std::uint32_t l, const std::uint32_t m) const
{
   double s = 2.0;

   // m = 0, no factor of two
   if (m == 0)
   {
      s = 1.0;
   }

   return s;
}

double ScaleLq::operator()(const std::uint32_t l) const
{
   double l_ = static_cast<double>(l);
   double s = std::pow(l_*(l_ + 1.0), 2);
   return s;
}

double ScaleLst::operator()(const std::uint32_t l) const
{
   double l_ = static_cast<double>(l);
   double s = l_*(l_ + 1.0);
   return s;
}

} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
