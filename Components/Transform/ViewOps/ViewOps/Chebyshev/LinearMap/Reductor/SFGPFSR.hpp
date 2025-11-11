/**
 * @file SFGPFR.hpp
 * @brief Chebyshev LinearMap reductor generic SpectralFftGridPointwiseFftSpectralReductor operator
 */
#pragma once

// External includes
//
#include <memory>

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
template <class Tout, class Tpower, class Tin, class SpecInBackend, class FftInBackend,
   class GridBackend, class PointBackend, class Functor, class FftOutBackend, class SpecOutBackend, class ReductorBackend>
class SFGPFSROp
    : public Operator::UnaryBaseOp<
         SFGPFSROp<Tout, Tpower, Tin, SpecInBackend, FftInBackend, GridBackend, PointBackend, Functor, FftOutBackend, SpecOutBackend, ReductorBackend>, Tout, Tin>
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
      SFGPFSROp<Tout, Tpower, Tin, SpecInBackend, FftInBackend, GridBackend, PointBackend, Functor, FftOutBackend, SpecOutBackend, ReductorBackend>, Tout, Tin>;
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
   class GridBackend, class PointBackend, class Functor, class FftOutBackend, class SpecOutBackend, class ReductorBackend>
SFGPFSROp<Tout, Tpower, Tin, SpecInBackend, FftInBackend, GridBackend, PointBackend, Functor, FftOutBackend, SpecOutBackend, ReductorBackend>::SFGPFSROp(
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
   class GridBackend, class PointBackend, class Functor, class FftOutBackend, class SpecOutBackend, class ReductorBackend>
void SFGPFSROp<Tout, Tpower, Tin, SpecInBackend, FftInBackend, GridBackend, PointBackend, Functor, FftOutBackend, SpecOutBackend, ReductorBackend>::applyImpl(
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
         _tmpData = std::move(Memory::MemBlock<typename Tin::ScalarType>(
            in.size(), _mem.get()));
         _tmpView = Tin(_tmpData.data(), _tmpData.size(), in.dims(),
            in.pointers(), in.indices(), in.lds());
      }
   }

   // setup power storage
   if (_powerView.data() == nullptr)
   {
      _powerData = std::move(Memory::MemBlock<typename Tpower::ScalarType>(
               in.size(), _mem.get()));
      _powerView = Tpower(_powerData.data(), _powerData.size(), in.dims(),
            in.pointers(), in.indices(), in.lds());
   }

   std::cerr << "IN:" << std::endl;
   for(int i = 0; i < in.size(); i++)
   {
      std::cerr << in.data()[i] << std::endl;
   }

   // spectral operation on input
   mSpecIn->apply(_tmpView, in, 1.0);

   std::cerr << "SPECIN OUT:" << std::endl;
   for(int i = 0; i < _tmpView.size(); i++)
   {
      std::cerr << _tmpView.data()[i] << std::endl;
   }

   // FFT
   mFftIn->apply(_tmpView, _tmpView);

   std::cerr << "FFT BWD OUT:" << std::endl;
   for(int i = 0; i < _tmpView.size(); i++)
   {
      std::cerr << _tmpView.data()[i] << std::endl;
   }

   // grid operation
   mGrid->apply(_tmpView, _tmpView, 1);

   std::cerr << "GRID OUT:" << std::endl;
   for(int i = 0; i < _tmpView.size(); i++)
   {
      std::cerr << _tmpView.data()[i] << std::endl;
   }

   // pointwise operation
   mPoint->apply(_powerView, _tmpView);

   std::cerr << "POINT OUT:" << std::endl;
   for(int i = 0; i < _powerView.size(); i++)
   {
      std::cerr << _powerView.data()[i] << std::endl;
   }

   // FFT
   mFftOut->apply(_powerView, _powerView);

   std::cerr << "FFT OUT:" << std::endl;
   for(int i = 0; i < _powerView.size(); i++)
   {
      std::cerr << _powerView.data()[i] << std::endl;
   }

   // spectral operation on input
   ScaleType fftScaling = 1.0 / static_cast<ScaleType>(2 * _powerView.dims()[0]);
   mSpecOut->apply(_powerView, _powerView, fftScaling);

   std::cerr << "SPECOUT OUT:" << std::endl;
   for(int i = 0; i < _powerView.size(); i++)
   {
      std::cerr << _powerView.data()[i] << std::endl;
   }

   // reduction operation
   mReductor->apply(out, _powerView);

   std::cerr << "REDUCTION OUT:" << std::endl;
   for(int i = 0; i < out.size(); i++)
   {
      std::cerr << out.data()[i] << std::endl;
   }
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
