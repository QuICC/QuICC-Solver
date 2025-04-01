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
#include "Fft/FftTags.hpp"
#include "ViewOps/Chebyshev/LinearMap/FftTypeMap.hpp"
#include "ViewOps/Chebyshev/LinearMap/Grid.hpp"
#include "ViewOps/Chebyshev/LinearMap/Integrator/GFS.hpp"
#include "ViewOps/Chebyshev/LinearMap/Projector/SFG.hpp"
#include "ViewOps/Chebyshev/LinearMap/Spec.hpp"
#include "ViewOps/Chebyshev/LinearMap/Tags.hpp"


namespace QuICC {
namespace Transform {
namespace Chebyshev {
namespace LinearMap {

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

template <class Backend, class Tmods, class Operation, std::uint16_t Treatment>
struct Spec;

template <class Backend, class Tmods, class Operation, std::uint16_t Treatment>
using Spec_t = typename Spec<Backend, Tmods, Operation, Treatment>::type;

template <class Tmods, class Operation, std::uint16_t Treatment>
struct Spec<viewCpu_t, Tmods, Operation, Treatment>
{
   using type = typename Cpu::SpecOp<Tmods, Tmods, Operation, Treatment>;
};

#ifdef QUICC_HAS_CUDA_BACKEND
template <class Tmods, class Operation, std::uint16_t Treatment>
struct Spec<viewGpu_t, Tmods, Operation, Treatment>
{
   using type = typename Cuda::SpecOp<Tmods, Tmods, Operation, Treatment>;
};
#endif

#ifdef QUICC_USE_VKFFT
template <class Tmods, class Operation, std::uint16_t Treatment>
struct Spec<viewGpuVkFFT_t, Tmods, Operation, Treatment>
{
   using type = typename Cuda::SpecOp<Tmods, Tmods, Operation, Treatment>;
};
#endif

template <class Backend, class Tphys, class Operation, std::uint16_t Treatment>
struct Grid;

template <class Backend, class Tphys, class Operation, std::uint16_t Treatment>
using Grid_t = typename Grid<Backend, Tphys, Operation, Treatment>::type;

template <class Tphys, class Operation, std::uint16_t Treatment>
struct Grid<viewCpu_t, Tphys, Operation, Treatment>
{
   using type = typename Cpu::GridOp<Tphys, Tphys, Operation, Treatment>;
};

#ifdef QUICC_HAS_CUDA_BACKEND
template <class Tphys, class Operation, std::uint16_t Treatment>
struct Grid<viewGpu_t, Tphys, Operation, Treatment>
{
   using type = typename Cuda::GridOp<Tphys, Tphys, Operation, Treatment>;
};
#endif

#ifdef QUICC_USE_VKFFT
template <class Tphys, class Operation, std::uint16_t Treatment>
struct Grid<viewGpuVkFFT_t, Tphys, Operation, Treatment>
{
   using type = typename Cuda::GridOp<Tphys, Tphys, Operation, Treatment>;
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
   using backendGrid_t = Grid_t<BACKEND, Tin, grid_id, none_t>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type2_t>;
   using backendSpec_t = Spec_t<BACKEND, Tout, spec_id, ndealias_out>;
   using type =
      Integrator::GFSOp<Tout, Tin, backendGrid_t, backendFft_t, backendSpec_t>;
   ;
};

/// @brief Op Y1 type map
/// Integrator of Y^1
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, Y1_t, fwd_t, BACKEND>
{
   using backendGrid_t = Grid_t<BACKEND, Tin, grid_id, none_t>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type2_t>;
   using backendSpec_t = Spec_t<BACKEND, Tout, spec_y1, ndealias_out>;
   using type =
      Integrator::GFSOp<Tout, Tin, backendGrid_t, backendFft_t, backendSpec_t>;
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
   using backendSpec_t = Spec_t<BACKEND, Tin, spec_id, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_id, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

/// @brief Op DivY1 type map
/// Projector of 1/Y
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, DivY1_t, bwd_t, BACKEND>
{
   using backendSpec_t = Spec_t<BACKEND, Tin, spec_id, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_divy1, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

/// @brief Op DivY2 type map
/// Projector of 1/Y^2
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, DivY2_t, bwd_t, BACKEND>
{
   using backendSpec_t = Spec_t<BACKEND, Tin, spec_id, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_divy2, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

/// @brief Op D1 type map
/// Projector of D^1
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D1_t, bwd_t, BACKEND>
{
   using backendSpec_t = Spec_t<BACKEND, Tin, spec_d1, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_id, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

/// @brief Op D2 type map
/// Projector of D^2
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D2_t, bwd_t, BACKEND>
{
   using backendSpec_t = Spec_t<BACKEND, Tin, spec_d2, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_id, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

/// @brief Op D3 type map
/// Projector of D^3
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D3_t, bwd_t, BACKEND>
{
   using backendSpec_t = Spec_t<BACKEND, Tin, spec_d3, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_id, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

/// @brief Op D4 type map
/// Projector of D^4
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D4_t, bwd_t, BACKEND>
{
   using backendSpec_t = Spec_t<BACKEND, Tin, spec_d4, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_id, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

/// @brief Op D1Y1 type map
/// Projector of D^1Y1
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, D1Y1_t, bwd_t, BACKEND>
{
   using backendSpec_t =
      Spec_t<BACKEND, Tin, spec_d1y1, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_id, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

/// @brief Op DivY1D1Y1 type map
/// Projector of 1/Y^1 D^1Y1
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, DivY1D1Y1_t, bwd_t, BACKEND>
{
   using backendSpec_t =
      Spec_t<BACKEND, Tin, spec_d1y1, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_divy1, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

/// @brief Op SphRadLapl type map
/// Projector of 1/Y^2 D^1Y2D^1
/// @tparam Tout
/// @tparam Tin
/// @tparam BACKEND
template <class Tout, class Tin, class BACKEND>
struct OpsTypeMap<Tout, Tin, SphRadLapl_t, bwd_t, BACKEND>
{
   using backendSpec_t =
      Spec_t<BACKEND, Tin, spec_d1y2d1, ndealias_in | zero_pad>;
   using backendFft_t =
      details::Fft_t<BACKEND, Tout, Tin, QuICC::Fft::dct_type3_t>;
   using backendGrid_t = Grid_t<BACKEND, Tout, grid_divy2, none_t>;
   using type =
      Projector::SFGOp<Tout, Tin, backendSpec_t, backendFft_t, backendGrid_t>;
   ;
};

} // namespace details

} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
