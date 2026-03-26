/**
 * @file GFS.hpp
 * @brief Chebyshev LinearMap integrator generic GridFftSpectral operator
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
//
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Operator/Binary.hpp"
#include "Operator/Unary.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Transform {
namespace Chebyshev {
namespace LinearMap {
/// @brief namespace for Chebyshev LinearMap integrator (physical to modal
/// space)
namespace Integrator {

/// @brief This class implements a Chebyshev LinearMap grid operation, Fft,
/// spectral operation integration
/// @tparam Tout output physical space type
/// @tparam Tin input modes type
/// @tparam GridBackend type of grid operator
/// @tparam FftBackend  type of FFT operator
/// @tparam SpecBackend type of spectral operator
template <class Tout, class Tin, class GridBackend, class FftBackend,
   class SpecBackend>
class GFSOp
    : public Operator::UnaryBaseOp<
         GFSOp<Tout, Tin, GridBackend, FftBackend, SpecBackend>, Tout, Tin>
{
public:
   /// @brief type of scale parameter, i.e. float 32/64 bits
   using ScaleType = double;
   /// @brief constructor with user defined scaling factor
   /// @param mem
   /// @param lower lower bound
   /// @param upper upper bound
   /// @param scale
   GFSOp(std::shared_ptr<Memory::memory_resource> mem, const double lower,
      const double upper, ScaleType scale = 1.0);
   /// @brief default constructor
   GFSOp() = delete;
   /// @brief dtor
   ~GFSOp() = default;

private:
   /// @brief action implementation that does not overwrite the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   void applyImpl(Tout& out, const Tin& in);
   /// @brief pointer to grid operator
   std::unique_ptr<Operator::BinaryOp<Tin, Tin, ScaleType>> mGrid;
   /// @brief pointer to FFT operator
   std::unique_ptr<Operator::UnaryOp<Tout, Tin>> mFft;
   /// @brief pointer to spectral operator
   std::unique_ptr<Operator::BinaryOp<Tin, Tin, ScaleType>> mSpec;
   /// @brief give access to base class
   friend Operator::UnaryBaseOp<
      GFSOp<Tout, Tin, GridBackend, FftBackend, SpecBackend>, Tout, Tin>;
   /// @brief memory resource
   /// needs shared ptr for memory pools
   /// note, this must call the dtor last
   /// otherwise we cannot dealloc data
   std::shared_ptr<Memory::memory_resource> _mem;
   /// @brief temporary memory block
   Memory::MemBlock<typename Tin::ScalarType> _tmpData;
   /// @brief Input View for the operator
   Tin _tmpView;
};

template <class Tout, class Tin, class GridBackend, class FftBackend,
   class SpecBackend>
GFSOp<Tout, Tin, GridBackend, FftBackend, SpecBackend>::GFSOp(
   std::shared_ptr<Memory::memory_resource> mem, const double lower,
   const double upper, ScaleType scale) :
    mGrid(std::make_unique<GridBackend>(lower, upper)),
    mFft(std::make_unique<FftBackend>()),
    mSpec(std::make_unique<SpecBackend>(lower, upper)),
    _mem(mem)
{}

template <class Tout, class Tin, class GridBackend, class FftBackend,
   class SpecBackend>
void GFSOp<Tout, Tin, GridBackend, FftBackend, SpecBackend>::applyImpl(
   Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix(
      "Chebyshev::LinearMap::Projector::GFSOp::applyImpl");

   // setup tmp storage
   if (_tmpView.data() == nullptr)
   {
      //        // special treatment for the pure integrator is not necessary
      //        // but it avoids an allocation and copy
      //        if constexpr (DiffBackend::TreatmentValue == none_m &&
      //        DiffBackend::OrderValue == 0)
      //        {
      //            _tmpView = in;
      //        }
      //        else
      {
         _tmpData = std::move(
            Memory::MemBlock<typename Tin::ScalarType>(in.size(), _mem.get()));
         _tmpView = Tin(_tmpData.data(), _tmpData.size(), in.dims(),
            in.pointers(), in.indices(), in.lds());
      }
   }

   // grid operation
   mGrid->apply(_tmpView, in, 1.0);

   // FFT
   mFft->apply(_tmpView, _tmpView);

   // spectral operation
   ScaleType fftScaling = 1.0 / static_cast<ScaleType>(2 * _tmpView.dims()[0]);
   mSpec->apply(out, _tmpView, fftScaling);
}

} // namespace Integrator
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
