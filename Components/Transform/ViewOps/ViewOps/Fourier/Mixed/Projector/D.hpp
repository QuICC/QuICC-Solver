/**
 * @file D.hpp
 * @brief Mixed projector diff operator
 */
#pragma once

// External includes
//
#include <memory>

// Project includes
//
#include "Memory/Memory.hpp"
#include "Memory/MemoryResource.hpp"
#include "Operator/Unary.hpp"
#include "Operator/Binary.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Mixed {
/// @brief namespace for mixed projectors (modal to physical space)
namespace Projector {

using namespace QuICC::Operator;
using namespace QuICC::Memory;

/// @brief This class implements a Fourier differentiation and projection
/// @tparam Tout output physical space type
/// @tparam Tin input modes type
/// @tparam FftBackend  type of FFT operator
/// @tparam DiffBackend type of operator
template<class Tout, class Tin, class FftBackend, class DiffBackend>
class DOp : public UnaryBaseOp<DOp<Tout, Tin, FftBackend, DiffBackend>, Tout, Tin> {
public:
    /// @brief type of scale parameter, i.e. float 32/64 bits
    using ScaleType = typename Tout::ScalarType;
    /// @brief constructor with user defined scaling factor
    /// @param mem
    /// @param scale
    DOp(std::shared_ptr<memory_resource> mem,
      ScaleType scale = 1.0);
    /// @brief default constructor
    DOp() = delete;
    /// @brief dtor
    ~DOp() = default;
private:
    /// @brief action implementation that does not overwrite the input
    /// @param out differentiatied physical space coefficient
    /// @param in input modes
    void applyImpl(Tout& out, const Tin& in);
    /// @brief pointer to FFT operator
    std::unique_ptr<UnaryOp<Tout, Tin>> mFft;
    /// @brief pointer to differentiation operator
    std::unique_ptr<BinaryOp<Tin, Tin, ScaleType>> mDiff;
    /// @brief give access to base class
    friend UnaryBaseOp<DOp<Tout, Tin, FftBackend, DiffBackend>, Tout, Tin>;
   /// @brief memory resource
   /// needs shared ptr for memory pools
   /// note, this must call the dtor last
   /// otherwise we cannot dealloc data
   std::shared_ptr<memory_resource> _mem;
   /// @brief temporary memory block
   MemBlock<typename Tin::ScalarType> _tmpData;
   /// @brief View for the operator
   Tin _tmpView;
};

template<class Tout, class Tin, class FftBackend, class DiffBackend>
DOp<Tout, Tin, FftBackend, DiffBackend>::DOp(std::shared_ptr<memory_resource> mem, ScaleType scale) : mFft(std::make_unique<FftBackend>()),
    mDiff(std::make_unique<DiffBackend>(scale)), _mem(mem)
{
}

template<class Tout, class Tin, class FftBackend, class DiffBackend>
void DOp<Tout, Tin, FftBackend, DiffBackend>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("Fourier::Mixed::Projector::DOp::applyImpl");

    // setup tmp storage
    if (_tmpView.data() == nullptr)
    {
        // special treatment for the pure projector is not necessary
        // but it avoids an allocation and copy
        if constexpr (DiffBackend::TreatmentValue == none_m && DiffBackend::OrderValue == 0)
        {
            _tmpView = in;
        }
        else
        {
        _tmpData = std::move(QuICC::Memory::MemBlock<typename Tin::ScalarType>(in.size(), _mem.get()));
        _tmpView = Tin(_tmpData.data(), _tmpData.size(), in.dims(), in.pointers(), in.indices(), in.lds());
        }
    }

    // differentiate
    mDiff->apply(_tmpView, in, 1.0);

    // FFT
    mFft->apply(out, _tmpView);
}

} // namespace Projector
} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
