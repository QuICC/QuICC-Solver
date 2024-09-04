/**
 * @file OpsTypes.hpp
 * @brief Mapping Generic Fourier operator types to specific Ops
 */

#pragma once

// External includes
//
#include <cstdint>

// Project includes
//
#include "Fft/Fft.hpp"
#include "ViewOps/Fourier/FftTypeMap.hpp"
#include "ViewOps/Fourier/Mixed/Diff.hpp"
#include "ViewOps/Fourier/Mixed/Integrator/D.hpp"
#include "ViewOps/Fourier/Mixed/Projector/D.hpp"
#include "ViewOps/Fourier/Tags.hpp"


namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Mixed {

/// @brief This namespace hides implementation details
namespace details {

/// @brief Generic mapping
/// @tparam Tout
/// @tparam Tin
/// @tparam TAG kind
/// @tparam DIR fwd_t or bwd_t
/// @tparam BACKEND
template <class Tout, class Tin, class TAG, class DIR, class BACKEND>
struct OpsTypeMap
{
   using type = void;
};

} // namespace details

/// @brief Convenience wrapper
/// (you cannot specialize type aliases)
/// @tparam Tout
/// @tparam Tin
/// @tparam TAG kind
/// @tparam DIR fwd_t or bwd_t
/// @tparam BACKEND
template <class Tout, class Tin, class TAG, class DIR, class BACKEND>
using OpsType =
   typename details::OpsTypeMap<Tout, Tin, TAG, DIR, BACKEND>::type;


namespace details {

template <class Backend, class Tmods, std::size_t Order, class Direction,
   std::uint16_t Treatment>
struct Diff;

template <class Backend, class Tmods, std::size_t Order, class Direction,
   std::uint16_t Treatment>
using Diff_t = typename Diff<Backend, Tmods, Order, Direction, Treatment>::type;

template <class Tmods, std::size_t Order, class Direction,
   std::uint16_t Treatment>
struct Diff<viewCpu_t, Tmods, Order, Direction, Treatment>
{
   using type = typename Cpu::DiffOp<Tmods, Tmods, Order, Direction, Treatment>;
};

#ifdef QUICC_HAS_CUDA_BACKEND
template <class Tmods, std::size_t Order, class Direction,
   std::uint16_t Treatment>
struct Diff<viewGpu_t, Tmods, Order, Direction, Treatment>
{
   using type =
      typename Cuda::DiffOp<Tmods, Tmods, Order, Direction, Treatment>;
};
#endif

#ifdef QUICC_USE_VKFFT
template <class Tmods, std::size_t Order, class Direction,
   std::uint16_t Treatment>
struct Diff<viewGpuVkFFT_t, Tmods, Order, Direction, Treatment>
{
   using type =
      typename Cuda::DiffOp<Tmods, Tmods, Order, Direction, Treatment>;
};
#endif

/// @brief Op P type map
/// Integrator only
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, P_t, fwd_t, BACKEND>
{
   using backendFft_t = Fourier::details::Fft_t<BACKEND, Tout, Tin>;
   using backendDiff_t =
      Diff_t<BACKEND, Tout, 0, fwd_t, QuICC::Transform::Fourier::none_m>;
   using type = Integrator::DOp<Tout, Tin, backendFft_t, backendDiff_t>;
   ;
};

/// @brief Op P type map
/// Projector only
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, P_t, bwd_t, BACKEND>
{
   using backendFft_t = Fourier::details::Fft_t<BACKEND, Tout, Tin>;
   using backendDiff_t =
      Diff_t<BACKEND, Tin, 0, bwd_t, QuICC::Transform::Fourier::none_m>;
   using type = Projector::DOp<Tout, Tin, backendFft_t, backendDiff_t>;
   ;
};

/// @brief Op D1 type map
/// Integrator only
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D1_t, fwd_t, BACKEND>
{
   using backendFft_t = Fourier::details::Fft_t<BACKEND, Tout, Tin>;
   using backendDiff_t =
      Diff_t<BACKEND, Tout, 1, fwd_t, QuICC::Transform::Fourier::none_m>;
   using type = Integrator::DOp<Tout, Tin, backendFft_t, backendDiff_t>;
   ;
};

/// @brief Op D1 type map
/// Projector only
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D1_t, bwd_t, BACKEND>
{
   using backendFft_t = Fourier::details::Fft_t<BACKEND, Tout, Tin>;
   using backendDiff_t =
      Diff_t<BACKEND, Tin, 1, bwd_t, QuICC::Transform::Fourier::none_m>;
   using type = Projector::DOp<Tout, Tin, backendFft_t, backendDiff_t>;
   ;
};

/// @brief Op D2 type map
/// Integrator only
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D2_t, fwd_t, BACKEND>
{
   using backendFft_t = Fourier::details::Fft_t<BACKEND, Tout, Tin>;
   using backendDiff_t =
      Diff_t<BACKEND, Tout, 2, fwd_t, QuICC::Transform::Fourier::none_m>;
   using type = Integrator::DOp<Tout, Tin, backendFft_t, backendDiff_t>;
   ;
};

/// @brief Op D2 type map
/// Projector only
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D2_t, bwd_t, BACKEND>
{
   using backendFft_t = Fourier::details::Fft_t<BACKEND, Tout, Tin>;
   using backendDiff_t =
      Diff_t<BACKEND, Tin, 2, bwd_t, QuICC::Transform::Fourier::none_m>;
   using type = Projector::DOp<Tout, Tin, backendFft_t, backendDiff_t>;
   ;
};

/// @brief Op D3 type map
/// Projector only
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D3_t, bwd_t, BACKEND>
{
   using backendFft_t = Fourier::details::Fft_t<BACKEND, Tout, Tin>;
   using backendDiff_t =
      Diff_t<BACKEND, Tin, 3, bwd_t, QuICC::Transform::Fourier::none_m>;
   using type = Projector::DOp<Tout, Tin, backendFft_t, backendDiff_t>;
   ;
};

/// @brief Op D1_P type map
/// Integrator only
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D1_P_t, fwd_t, BACKEND>
{
   using backendFft_t = Fourier::details::Fft_t<BACKEND, Tout, Tin>;
   using backendDiff_t =
      Diff_t<BACKEND, Tout, 1, fwd_t, QuICC::Transform::Fourier::zeroP_m>;
   using type = Integrator::DOp<Tout, Tin, backendFft_t, backendDiff_t>;
   ;
};

/// @brief Op D1_Neg type map
/// Integrator only
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D1_Neg_t, fwd_t, BACKEND>
{
   using backendFft_t = Fourier::details::Fft_t<BACKEND, Tout, Tin>;
   using backendDiff_t =
      Diff_t<BACKEND, Tout, 1, fwd_t, QuICC::Transform::Fourier::zeroMinusP_m>;
   using type = Integrator::DOp<Tout, Tin, backendFft_t, backendDiff_t>;
   ;
};

} // namespace details

} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
