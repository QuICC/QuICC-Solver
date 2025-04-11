/**
 * @file SFGFS.hpp
 * @brief Chebyshev LinearMap reductor generic SpectralFftGridFftSpectral operator
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
/// @brief namespace for Chebyshev LinearMap reductor (modal to energy)
namespace Reductor {

/// @brief This class implements a Chebyshev LinearMap Spectral space operation,
/// FFT and grid operation
/// @tparam Tout output energy stype
/// @tparam Tin input modes type
/// @tparam SpecInBackend type of spectral operator
/// @tparam FftInBackend  type of FFT operator
/// @tparam GridBackend type of grid operator
/// @tparam FftOutBackend  type of FFT operator
/// @tparam SpecOutBackend type of spectral operator
template <class Tout, class Tin, class SpecInBackend, class FftInBackend,
   class GridBackend, class FftOutBackend, class SpecOutBackend>
class SFGFSOp
    : public Operator::UnaryBaseOp<
         SFGFSOp<Tout, Tin, SpecInBackend, FftInBackend, GridBackend, FftOutBackend, SpecOutBackend>, Tout, Tin>
{
public:
   /// @brief type of scale parameter, i.e. float 32/64 bits
   using ScaleType = double;
   /// @brief constructor with user defined scaling factor
   /// @param mem
   /// @param lower Upper bound
   /// @param upper lower bound
   /// @param scale
   SFGFSOp(std::shared_ptr<Memory::memory_resource> mem, const double lower,
      const double upper, ScaleType scale = 1.0);
   /// @brief default constructor
   SFGFSOp() = delete;
   /// @brief dtor
   ~SFGFSOp() = default;

private:
   /// @brief action implementation that does not overwrite the input
   /// @param out differentiatied physical space coefficient
   /// @param in input modes
   void applyImpl(Tout& out, const Tin& in);
   /// @brief pointer to spectral operator
   std::unique_ptr<Operator::BinaryOp<Tin, Tin, ScaleType>> mSpecIn;
   /// @brief pointer to FFT operator
   std::unique_ptr<Operator::UnaryOp<Tout, Tin>> mFftIn;
   /// @brief pointer to grid operator
   std::unique_ptr<Operator::BinaryOp<Tout, Tout, ScaleType>> mGrid;
   /// @brief pointer to FFT operator
   std::unique_ptr<Operator::UnaryOp<Tout, Tin>> mFftOut;
   /// @brief pointer to spectral operator
   std::unique_ptr<Operator::BinaryOp<Tin, Tin, ScaleType>> mSpecOut;
   /// @brief give access to base class
   friend Operator::UnaryBaseOp<
      SFGFSOp<Tout, Tin, SpecInBackend, FftInBackend, GridBackend, FftOutBackend, SpecOutBackend>, Tout, Tin>;
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

template <class Tout, class Tin, class SpecInBackend, class FftInBackend,
   class GridBackend, class FftOutBackend, class SpecOutBackend>
SFGFSOp<Tout, Tin, SpecInBackend, FftInBackend, GridBackend, FftOutBackend, SpecOutBackend>::SFGFSOp(
   std::shared_ptr<Memory::memory_resource> mem, const double lower,
   const double upper, ScaleType scale) :
    mSpecIn(std::make_unique<SpecInBackend>(lower, upper)),
    mFftIn(std::make_unique<FftInBackend>()),
    mGrid(std::make_unique<GridBackend>(lower, upper)),
    mFftOut(std::make_unique<FftInBackend>()),
    mSpecOut(std::make_unique<SpecInBackend>(lower, upper)),
    _mem(mem)
{}

template <class Tout, class Tin, class SpecInBackend, class FftInBackend,
   class GridBackend, class FftOutBackend, class SpecOutBackend>
void SFGFSOp<Tout, Tin, SpecInBackend, FftInBackend, GridBackend, FftOutBackend, SpecOutBackend>::applyImpl(
   Tout& out, const Tin& in)
{
   Profiler::RegionFixture<4> fix(
      "Chebyshev::LinearMap::Projector::SFGFSOp::applyImpl");

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
         _tmpData = std::move(Memory::MemBlock<typename Tin::ScalarType>(
            in.size(), _mem.get()));
         _tmpView = Tout(_tmpData.data(), _tmpData.size(), out.dims(),
            out.pointers(), out.indices(), out.lds());
      }
   }

   // spectral operation
   mSpecIn->apply(_tmpView, in, 1.0);

   // FFT
   mFftIn->apply(_tmpView, _tmpView);

   // grid operation
   mGrid->apply(out, _tmpView, 1.0);

   // FFT
   mFftOut->apply(_tmpView, _tmpView);

   // spectral operation
   mSpecOut->apply(_tmpView, in, 1.0);
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
