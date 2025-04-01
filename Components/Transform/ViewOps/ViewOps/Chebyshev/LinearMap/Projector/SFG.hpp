/**
 * @file SFG.hpp
 * @brief Chebyshev LinearMap projector generic SpectralFftGrid operator
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
/// @brief namespace for Chebyshev LinearMap projectors (modal to physical
/// space)
namespace Projector {

/// @brief This class implements a Chebyshev LinearMap Spectral space operation,
/// FFT projection and grid operation
/// @tparam Tout output physical space type
/// @tparam Tin input modes type
/// @tparam SpecBackend type of spectral operator
/// @tparam FftBackend  type of FFT operator
/// @tparam GridBackend type of grid operator
template <class Tout, class Tin, class SpecBackend, class FftBackend,
   class GridBackend>
class SFGOp
    : public Operator::UnaryBaseOp<
         SFGOp<Tout, Tin, SpecBackend, FftBackend, GridBackend>, Tout, Tin>
{
public:
   /// @brief type of scale parameter, i.e. float 32/64 bits
   using ScaleType = double;
   /// @brief constructor with user defined scaling factor
   /// @param mem
   /// @param lower Upper bound
   /// @param upper lower bound
   /// @param scale
   SFGOp(std::shared_ptr<Memory::memory_resource> mem, const double lower,
      const double upper, ScaleType scale = 1.0);
   /// @brief default constructor
   SFGOp() = delete;
   /// @brief dtor
   ~SFGOp() = default;

private:
   /// @brief action implementation that does not overwrite the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   void applyImpl(Tout& out, const Tin& in);
   /// @brief pointer to spectral operator
   std::unique_ptr<Operator::BinaryOp<Tin, Tin, ScaleType>> mSpec;
   /// @brief pointer to FFT operator
   std::unique_ptr<Operator::UnaryOp<Tout, Tin>> mFft;
   /// @brief pointer to grid operator
   std::unique_ptr<Operator::BinaryOp<Tout, Tout, ScaleType>> mGrid;
   /// @brief give access to base class
   friend Operator::UnaryBaseOp<
      SFGOp<Tout, Tin, SpecBackend, FftBackend, GridBackend>, Tout, Tin>;
   /// @brief memory resource
   /// needs shared ptr for memory pools
   /// note, this must call the dtor last
   /// otherwise we cannot dealloc data
   std::shared_ptr<Memory::memory_resource> _mem;
   /// @brief temporary memory block
   Memory::MemBlock<typename Tin::ScalarType> _tmpData;
   /// @brief View for the operator
   Tin _tmpView;
};

template <class Tout, class Tin, class SpecBackend, class FftBackend,
   class GridBackend>
SFGOp<Tout, Tin, SpecBackend, FftBackend, GridBackend>::SFGOp(
   std::shared_ptr<Memory::memory_resource> mem, const double lower,
   const double upper, ScaleType scale) :
    mSpec(std::make_unique<SpecBackend>(lower, upper)),
    mFft(std::make_unique<FftBackend>()),
    mGrid(std::make_unique<GridBackend>(lower, upper)),
    _mem(mem)
{}

template <class Tout, class Tin, class SpecBackend, class FftBackend,
   class GridBackend>
void SFGOp<Tout, Tin, SpecBackend, FftBackend, GridBackend>::applyImpl(
   Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix(
      "Chebyshev::LinearMap::Projector::SFGOp::applyImpl");

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
         _tmpData = std::move(Memory::MemBlock<typename Tout::ScalarType>(
            out.size(), _mem.get()));
         _tmpView = Tout(_tmpData.data(), _tmpData.size(), out.dims(),
            out.pointers(), out.indices(), out.lds());
      }
   }

   // spectral operation
   mSpec->apply(_tmpView, in, 1.0);

   // FFT
   mFft->apply(_tmpView, _tmpView);

   // grid operation
   mGrid->apply(out, _tmpView, 1.0);
}

} // namespace Projector
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
