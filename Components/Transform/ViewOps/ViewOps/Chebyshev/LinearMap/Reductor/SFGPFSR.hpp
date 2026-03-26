/**
 * @file SFGPFR.hpp
 * @brief Chebyshev LinearMap reductor generic SpectralFftGridPointwiseFftSpectralReductor operator
 */
#pragma once

// External includes
//
#include <memory>
#include <cstdint>

// Project includes
//
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Operator/Nary.hpp"
#include "Operator/Binary.hpp"
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"

#include <iostream>
namespace QuICC {
namespace Transform {
namespace Chebyshev {
namespace LinearMap {
/// @brief namespace for Chebyshev LinearMap reductor (modal to energy)
namespace Reductor {

/// @brief This class implements a Chebyshev LinearMap Spectral space operation,
/// FFT and grid operation
/// @tparam Tout output energy type
/// @tparam Tpower power type
/// @tparam Tin input modes type
/// @tparam SpecInBackend type of spectral operator on input
/// @tparam FftInBackend  type of FFT operator on input
/// @tparam GridBackend type of grid operator
/// @tparam PointBackend type of pointwise operator
/// @tparam FftOutBackend  type of FFT operator output
/// @tparam SpecOutBackend type of spectral operator on output
/// @tparam ReductorBackend type of reductor operator
/// @tparam SizeMult multiplier for size
/// @tparam SizeAdd additional modes for size
template <class Tout, class Tpower, class Tin, class SpecInBackend, class FftInBackend,
   class GridBackend, class PointBackend, class Functor, class FftOutBackend, class SpecOutBackend, class ReductorBackend, std::uint32_t SizeMult, std::uint32_t SizeAdd>
class SFGPFSROp
    : public Operator::UnaryBaseOp<
         SFGPFSROp<Tout, Tpower, Tin, SpecInBackend, FftInBackend, GridBackend, PointBackend, Functor, FftOutBackend, SpecOutBackend, ReductorBackend, SizeMult, SizeAdd>, Tout, Tin>
{
public:
   /// @brief type of scale parameter, i.e. float 32/64 bits
   using ScaleType = double;
   /// @brief constructor with user defined scaling factor
   /// @param mem
   /// @param lower Upper bound
   /// @param upper lower bound
   /// @param scale
   SFGPFSROp(std::shared_ptr<Memory::memory_resource> mem, const double lower,
      const double upper, ScaleType scale = 1.0);
   /// @brief default constructor
   SFGPFSROp() = delete;
   /// @brief dtor
   ~SFGPFSROp() = default;

private:
   /// @brief action implementation that does not overwrite the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   void applyImpl(Tout& out, const Tin& in);
   /// @brief pointer to spectral operator
   std::unique_ptr<Operator::BinaryOp<Tin, Tin, ScaleType>> mSpecIn;
   /// @brief pointer to FFT operator
   std::unique_ptr<Operator::UnaryOp<Tin, Tin>> mFftIn;
   /// @brief pointer to grid operator
   std::unique_ptr<Operator::BinaryOp<Tin, Tin, ScaleType>> mGrid;
   /// @brief pointer to pointwise operator
   std::unique_ptr<Operator::NaryOp<Tpower, Tin>> mPoint;
   /// @brief pointer to FFT operator
   std::unique_ptr<Operator::UnaryOp<Tpower, Tpower>> mFftOut;
   /// @brief pointer to spectral operator
   std::unique_ptr<Operator::BinaryOp<Tpower, Tpower, ScaleType>> mSpecOut;
   /// @brief pointer to spectral operator
   std::unique_ptr<Operator::UnaryOp<Tout, Tpower>> mReductor;
   /// @brief give access to base class
   friend Operator::UnaryBaseOp<
      SFGPFSROp<Tout, Tpower, Tin, SpecInBackend, FftInBackend, GridBackend, PointBackend, Functor, FftOutBackend, SpecOutBackend, ReductorBackend, SizeMult, SizeAdd>, Tout, Tin>;
   /// @brief memory resource
   /// needs shared ptr for memory pools
   /// note, this must call the dtor last
   /// otherwise we cannot dealloc data
   std::shared_ptr<Memory::memory_resource> _mem;
   /// @brief temporary memory block
   Memory::MemBlock<typename Tpower::ScalarType> _powerData;
   /// @brief View for the operator
   Tpower _powerView;
   /// @brief temporary memory block
   Memory::MemBlock<typename Tin::ScalarType> _tmpData;
   /// @brief View for the operator
   Tin _tmpView;
};

template <class Tout, class Tpower, class Tin, class SpecInBackend, class FftInBackend,
   class GridBackend, class PointBackend, class Functor, class FftOutBackend, class SpecOutBackend, class ReductorBackend, std::uint32_t SizeMult, std::uint32_t SizeAdd>
SFGPFSROp<Tout, Tpower, Tin, SpecInBackend, FftInBackend, GridBackend, PointBackend, Functor, FftOutBackend, SpecOutBackend, ReductorBackend, SizeMult, SizeAdd>::SFGPFSROp(
   std::shared_ptr<Memory::memory_resource> mem, const double lower,
   const double upper, ScaleType scale) :
    mSpecIn(std::make_unique<SpecInBackend>(lower, upper)),
    mFftIn(std::make_unique<FftInBackend>()),
    mGrid(std::make_unique<GridBackend>(lower, upper)),
    mPoint(std::make_unique<PointBackend>(Functor())),
    mFftOut(std::make_unique<FftOutBackend>()),
    mSpecOut(std::make_unique<SpecOutBackend>(lower, upper)),
    mReductor(std::make_unique<ReductorBackend>()),
    _mem(mem)
{}

template <class Tout, class Tpower, class Tin, class SpecInBackend, class FftInBackend,
   class GridBackend, class PointBackend, class Functor, class FftOutBackend, class SpecOutBackend, class ReductorBackend, std::uint32_t SizeMult, std::uint32_t SizeAdd>
void SFGPFSROp<Tout, Tpower, Tin, SpecInBackend, FftInBackend, GridBackend, PointBackend, Functor, FftOutBackend, SpecOutBackend, ReductorBackend, SizeMult, SizeAdd>::applyImpl(
   Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix(
      "Chebyshev::LinearMap::Energy::SFGFSOp::applyImpl");

   // setup tmp storage
   if (_tmpView.data() == nullptr)
   {
      //        // special treatment for the pure projector is not necessary
      //        // but it avoids an allocation and copy
      //        if constexpr (DiffBackend::TreatmentValue == none_m &&
      //        DiffBackend::OrderValue == 0)
      //        {
      //            _tmpView = in;
      //        }
      //        else
      {
         auto sze = SizeMult*in.dims()[0] + SizeAdd;
         auto memsze = sze*in.pointers()[1][in.pointers()[1].size()-1];

         _tmpData = std::move(Memory::MemBlock<typename Tin::ScalarType>(
            memsze, _mem.get()));
         std::array<typename Tin::IndexType, 3> tmpDims;
         std::copy(in.dims(), in.dims() + 3, tmpDims.data());
         tmpDims[0] = sze;
         _tmpView = Tin(_tmpData.data(), _tmpData.size(), tmpDims.data(),
            in.pointers(), in.indices(), tmpDims[0]);
      }
   }

   // setup power storage
   if (_powerView.data() == nullptr)
   {
      auto sze = SizeMult*in.dims()[0] + SizeAdd;
      auto memsze = sze*in.pointers()[1][in.pointers()[1].size()-1];

      _powerData = std::move(Memory::MemBlock<typename Tpower::ScalarType>(
               memsze, _mem.get()));
      std::array<typename Tin::IndexType, 3> powerDims;
      std::copy(in.dims(), in.dims() + 3, powerDims.data());
      powerDims[0] = sze;
      _powerView = Tpower(_powerData.data(), _powerData.size(), powerDims.data(),
            in.pointers(), in.indices(), powerDims[0]);
   }

   // spectral operation on input
   mSpecIn->apply(_tmpView, in, 1.0);

   // FFT
   mFftIn->apply(_tmpView, _tmpView);

   // grid operation
   mGrid->apply(_tmpView, _tmpView, 1);

   // pointwise operation
   mPoint->apply(_powerView, _tmpView);

   // FFT
   mFftOut->apply(_powerView, _powerView);

   // spectral operation on input
   ScaleType fftScaling = 1.0 / static_cast<ScaleType>(2 * _powerView.dims()[0]);
   mSpecOut->apply(_powerView, _powerView, fftScaling);

   // reduction operation
   mReductor->apply(out, _powerView);
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
