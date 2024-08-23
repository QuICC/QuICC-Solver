/**
 * @file BackendsMap.hpp
 * @brief
 */
#pragma once

// External includes
//

// Project includes
//
#include "Graph/Tags.hpp"
#include "Fft/Fft.hpp"
#include "ViewOps/Fourier/Mixed/Projector/D.hpp"
#include "ViewOps/Fourier/Mixed/Integrator/D.hpp"
#include "ViewOps/Fourier/Mixed/Diff.hpp"
#include "ViewOps/Worland/Tags.hpp"
#include "ViewOps/Quadrature/Impl.hpp"
#include "ViewOps/Quadrature/Op.hpp"
#include "ViewOps/Pointwise/Pointwise.hpp"
#include "ViewOps/Pointwise/Functors.hpp"
#include "ViewOps/Transpose/Op.hpp"
#include "Fft/Fft.hpp"
#include "ViewOps/Fourier/Mixed/Diff.hpp"

namespace QuICC {
namespace Graph {

// ///
// /// Pointwise
// ///

// template <class Backend, class Functor, class Tout, class... Targs>
// struct Pointwise;

// template <class Backend, class Functor, class Tout, class... Targs>
// using Pointwise_t = typename Pointwise<Backend, Functor, Tout, Targs...>::type;

// template <class Functor, class Tout, class... Targs>
// struct Pointwise<viewCpu_t, Functor, Tout, Targs...>
// {
//     using type = typename QuICC::Pointwise::Cpu::Op<Functor, Tout, Targs...>;
// };

// #ifdef QUICC_HAS_CUDA_BACKEND
// template <class Functor, class Tout, class... Targs>
// struct Pointwise<viewGpu_t, Functor, Tout, Targs...>
// {
//     using type = typename QuICC::Pointwise::Cuda::Op<Functor, Tout, Targs...>;
// };
// #endif

///
/// Fourier
///

template <class Backend, class Tout, class Tin>
struct Fft;

template <class Backend, class Tout, class Tin>
using Fft_t = typename Fft<Backend, Tout, Tin>::type;

template <class Tout, class Tin>
struct Fft<viewCpu_t, Tout, Tin>
{
    using type = typename QuICC::Fft::Fftw::FftOp<Tout, Tin>;
};

#ifdef QUICC_HAS_CUDA_BACKEND
template <class Tout, class Tin>
struct Fft<viewGpu_t, Tout, Tin>
{
    using type = typename QuICC::Fft::CuFft::FftOp<Tout, Tin>;
};
#endif

template <class Backend, class Tmods, std::size_t Order, class Direction, std::uint16_t Treatment>
struct MixedDiff;

template <class Backend, class Tmods, std::size_t Order, class Direction, std::uint16_t Treatment>
using MixedDiff_t = typename MixedDiff<Backend, Tmods, Order, Direction, Treatment>::type;

template <class Tmods, std::size_t Order, class Direction, std::uint16_t Treatment>
struct MixedDiff<viewCpu_t, Tmods, Order, Direction, Treatment>
{
    using type = typename QuICC::Transform::Fourier::Mixed::Cpu::DiffOp<Tmods, Tmods, Order, Direction, Treatment>;
};

#ifdef QUICC_HAS_CUDA_BACKEND
template <class Tmods, std::size_t Order, class Direction, std::uint16_t Treatment>
struct MixedDiff<viewGpu_t, Tmods, Order, Direction, Treatment>
{
    using type = typename QuICC::Transform::Fourier::Mixed::Cuda::DiffOp<Tmods, Tmods, Order, Direction, Treatment>;
};
#endif

} // namespace Transform
} // namespace QuICC
