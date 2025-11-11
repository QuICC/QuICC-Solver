/**
 * @file Tags.hpp
 * @brief Tag types
 */
#pragma once

// System includes
//
#include <array>
#include <cstdint>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y1D1Y1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2D1Y1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4D1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y3D1Y1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I4Y3.hpp"

namespace QuICC {
namespace Transform {
namespace Chebyshev {
namespace LinearMap {

//
// Tags
//

/// @brief tag type for projection direction.
/// Forwards i.e. physical to modal (integrator)
struct fwd_t
{
};

/// @brief tag type for projection direction.
/// Backwards i.e. modal to physical (projector)
struct bwd_t
{
};

/// @brief No sparse operator op type tag
struct spec_no_sparseop
{
   typedef void SparseOpType;
};

/// @brief tag type for identity spectral operation
struct spec_id: public spec_no_sparseop
{
   constexpr static std::array<std::size_t, 1> p = {0};
   constexpr static std::array<std::size_t, 1> t = {0};
};

/// @brief tag type for multiplication by y^1
struct spec_y1: public spec_no_sparseop
{
   constexpr static std::array<std::size_t, 1> p = {0};
   constexpr static std::array<std::size_t, 1> t = {1};
};

/// @brief tag type for spectral derivative
struct spec_d1: public spec_no_sparseop
{
   constexpr static std::array<std::size_t, 1> p = {1};
   constexpr static std::array<std::size_t, 1> t = {0};
};

/// @brief tag type for spectral second derivative
struct spec_d2: public spec_no_sparseop
{
   constexpr static std::array<std::size_t, 1> p = {2};
   constexpr static std::array<std::size_t, 1> t = {0};
};

/// @brief tag type for spectral third derivative
struct spec_d3: public spec_no_sparseop
{
   constexpr static std::array<std::size_t, 1> p = {3};
   constexpr static std::array<std::size_t, 1> t = {0};
};

/// @brief tag type for spectral fourth derivative
struct spec_d4: public spec_no_sparseop
{
   constexpr static std::array<std::size_t, 1> p = {4};
   constexpr static std::array<std::size_t, 1> t = {0};
};

/// @brief tag type for spectral operator D^1Y^1
struct spec_d1y1: public spec_no_sparseop
{
   constexpr static std::array<std::size_t, 1> p = {1};
   constexpr static std::array<std::size_t, 1> t = {1};
};

/// @brief tag type for spectral operator D^1Y^2D^1
struct spec_d1y2d1: public spec_no_sparseop
{
   constexpr static std::array<std::size_t, 2> p = {1, 1};
   constexpr static std::array<std::size_t, 2> t = {0, 2};
};

/// @brief tag type for weighting with integral
struct spec_int: public spec_no_sparseop
{
   constexpr static std::array<std::size_t, 1> p = {0};
   constexpr static std::array<std::size_t, 1> t = {0};
};

/// @brief sparse operator op type tag
template <typename Type>
struct spec_sparseop
{
   typedef Type SparseOpType;

   constexpr static std::array<std::size_t, 1> p = {0};
   constexpr static std::array<std::size_t, 1> t = {0};
};

struct spec_i2: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I2>
{
   typedef ::QuICC::SparseSM::Chebyshev::LinearMap::I2D1 SparseMeanOpType;
};

struct spec_i2d1: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I2D1>
{
   typedef ::QuICC::SparseSM::Chebyshev::LinearMap::I2 SparseMeanOpType;
};

struct spec_i2y1d1y1: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I2Y1D1Y1>
{
};

struct spec_i2y1: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I2Y1>
{
};

struct spec_i2y2d1y1: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I2Y2D1Y1>
{
};

struct spec_i2y2: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I2Y2>
{
};

struct spec_i4: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I4>
{
};

struct spec_i4d1: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I4D1>
{
   typedef ::QuICC::SparseSM::Chebyshev::LinearMap::I2 SparseMeanOpType;
};

struct spec_i4y3d1y1: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I4Y3D1Y1>
{
};

struct spec_i4y3: public spec_sparseop<::QuICC::SparseSM::Chebyshev::LinearMap::I4Y3>
{
};

/// @brief tag type for identity grid operation
struct grid_id
{
};

/// @brief tag type for division by Y^1 grid operation
struct grid_divy1
{
};

/// @brief tag type for division by Y^1 grid operation
struct grid_divy2
{
};

/// view cpu implementation tag
struct viewCpu_t
{
};

/// view gpu implementation tag
struct viewGpu_t
{
};

/// view gpu VkFFT implementation tag
struct viewGpuVkFFT_t
{
};

/// @brief no special treatment
constexpr std::uint16_t none_t = 0;

/// @brief dealias size from input
constexpr std::uint16_t ndealias_in = 1 << 1;

/// @brief dealias size from output
constexpr std::uint16_t ndealias_out = 1 << 2;

/// @brief zero pad
constexpr std::uint16_t zero_pad = 1 << 3;

/// @brief zero l = 0 modes
constexpr std::uint16_t zero_l0 = 1 << 4;

/// @brief special mean op modes
constexpr std::uint16_t mean_op = 1 << 5;

/// @brief P op type tag
struct P_t
{
};

/// @brief P_Zero op type tag
struct P_Zero_t
{
};

/// @brief P op type tag
struct Y1_t
{
};

/// @brief Y1_Zero op type tag
struct Y1_Zero_t
{
};

/// @brief DivY1_t op type tag
struct DivY1_t
{
};

/// @brief DivY2_t op type tag
struct DivY2_t
{
};

/// @brief D1_t op type tag
struct D1_t
{
};

/// @brief D2_t op type tag
struct D2_t
{
};

/// @brief D3_t op type tag
struct D3_t
{
};

/// @brief D4_t op type tag
struct D4_t
{
};

/// @brief D1Y1_t op type tag
struct D1Y1_t
{
};

/// @brief DivY1D1Y1_t op type tag
struct DivY1D1Y1_t
{
};

/// @brief SphRadLapl_t op type tag
struct SphRadLapl_t
{
};

/// @brief I2 op type tag
struct I2_t
{
};

/// @brief I2D1 op type tag
struct I2D1_t
{
};

/// @brief I2D1_I2 op type tag
struct I2D1_I2_t
{
};

/// @brief I2Y1D1Y1_Zero op type tag
struct I2Y1D1Y1_Zero_t
{
};

/// @brief I2Y1_Zero op type tag
struct I2Y1_Zero_t
{
};

/// @brief I2Y2D1Y1_Zero op type tag
struct I2Y2D1Y1_Zero_t
{
};

/// @brief I2Y2_Zero op type tag
struct I2Y2_Zero_t
{
};

/// @brief I2_I2D1 op type tag
struct I2_I2D1_t
{
};

/// @brief I4 op type tag
struct I4_t
{
};

/// @brief I4D1 op type tag
struct I4D1_t
{
};

/// @brief I4D1_I2 op type tag
struct I4D1_I2_t
{
};

/// @brief I4Y3D1Y1_Zero op type tag
struct I4Y3D1Y1_Zero_t
{
};

/// @brief I4Y3_Zero op type tag
struct I4Y3_Zero_t
{
};

/// @brief Energy op type tag
struct Energy_t
{
};

/// @brief EnergyD1 op type tag
struct EnergyD1_t
{
};

/// @brief EnergyD1Y1 op type tag
struct EnergyD1Y1_t
{
};

/// @brief Energy op type tag
struct EnergyY2_t
{
};

/// @brief EnergySLaplR2 op type tag
struct EnergySLaplR2_t
{
};


} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
