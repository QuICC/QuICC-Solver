/**
 * @file D.hpp
 * @brief Complex projector diff operator
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
#include "Profiler/Interface.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
/// @brief namespace for complex operators
namespace Complex {
/// @brief namespace for complex projectors (modal to physical space)
namespace Projector {


/// @brief This class implements a Fourier differentiation and projection
/// @tparam Tout output physical space type
/// @tparam Tin input modes type
/// @tparam FftBackend  type of FFT operator
/// @tparam DiffBackend type of operator
template<class Tout, class Tin, class FftBackend, class DiffBackend>
class DOp : public Operator::UnaryBaseOp<DOp<Tout, Tin, FftBackend, DiffBackend>, Tout, Tin> {
public:
    /// @brief Type of scale parameter, i.e. float 32/64 bits
    using ScaleType = typename Tout::ScalarType::value_type;
    /// @brief constructor with user defined scaling factor
    /// @param mem
    /// @param scale
    DOp(std::shared_ptr<Memory::memory_resource> mem,
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
    std::unique_ptr<Operator::UnaryOp<Tout, Tin>> mFft;
    /// @brief pointer to differentiation operator
    std::unique_ptr<Operator::UnaryOp<Tin, Tin>> mDiff;
    /// @brief Give access to base class
    friend Operator::UnaryBaseOp<DOp<Tout, Tin, FftBackend, DiffBackend>, Tout, Tin>;
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

template<class Tout, class Tin, class FftBackend, class DiffBackend>
DOp<Tout, Tin, FftBackend, DiffBackend>::DOp(std::shared_ptr<Memory::memory_resource> mem, ScaleType scale) : mFft(std::make_unique<FftBackend>()),
    mDiff(std::make_unique<DiffBackend>(scale)), _mem(mem)
{
}

template<class Tout, class Tin, class FftBackend, class DiffBackend>
void DOp<Tout, Tin, FftBackend, DiffBackend>::applyImpl(Tout& out, const Tin& in)
{
    Profiler::RegionFixture<4> fix("Fourier::Complex::Projector::DOp::applyImpl");

    // setup tmp storage
    if (_tmpView.data() == nullptr)
    {
        _tmpData = std::move(Memory::MemBlock<typename Tin::ScalarType>(in.size(), _mem.get()));
        _tmpView = Tin(_tmpData.data(), _tmpData.size(), in.dims(), in.pointers(), in.indices(), in.lds());
    }

    // differentiate
    mDiff->apply(_tmpView, in);

    // FFT
    mFft->apply(out, _tmpView);
}

} // namespace Projector
} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
