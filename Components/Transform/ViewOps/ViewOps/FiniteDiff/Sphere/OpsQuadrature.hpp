/**
 * @file OpsQuadrature.hpp
 * @brief Mapping Generic Finite differences operator builders to specific Ops
 */

#pragma once

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "FiniteDiff/Sphere/UniformRadialGrid.hpp"
#include "ViewOps/FiniteDiff/Sphere/Tags.hpp"

namespace QuICC {
namespace Transform {
namespace FiniteDiff {
namespace Sphere {

/// @brief This namespace hides implementation details
namespace details {

/// @brief Default radial grid for sphere
/// @tparam TAG kind
template <class TAG> struct OpsQuadratureMap
{
   using type = ::QuICC::FiniteDiff::Sphere::UniformRadialGrid;
};

} // namespace details

/// @brief Convenience wrapper
/// (you cannot specialize type aliases)
/// @tparam TAG kind
template <class TAG>
using OpsQuadrature = typename details::OpsQuadratureMap<TAG>::type;

// Specialized mapping for reductors
namespace details {

using namespace QuICC::FiniteDiff::Sphere;

/// @brief Compute quadrature for grid needed for power
/// computations
/// @tparam SHIFT
/// @param igrid
/// @param iweights
/// @param gSize
template <int SHIFT>
void computePowerQuadrature(Internal::Array& igrid, const int gSize)
{
   int nrgSize = gSize + 2 * SHIFT;
   FiniteDiff::Sphere::UniformRadialGrid quad;
   quad.computeGrid(igrid, nrgSize);
}

/// @brief Generic helper in order to avoid having to pass
/// extra parameters to energy integrators
template <class OP> struct EnergyHelperQuadMap;

/// @brief Energy helper in order to avoid having to pass
/// extra parameters to energy integrators
template <> struct EnergyHelperQuadMap<Energy_t>
{
   void computeQuadrature(Internal::Array& igrid,
      const int gSize)
   {
      computePowerQuadrature<1>(igrid, gSize);
   }
};

/// @brief EnergyD1R1 helper in order to avoid having to pass
/// extra parameters to energy integrators
template <> struct EnergyHelperQuadMap<EnergyD1R1_t>
{
   void computeQuadrature(Internal::Array& igrid, const int gSize)
   {
      computePowerQuadrature<0>(igrid, gSize);
   }
};

/// @brief EnergyR2 helper in order to avoid having to pass
/// extra parameters to energy integrators
template <> struct EnergyHelperQuadMap<EnergyR2_t>
{
   void computeQuadrature(Internal::Array& igrid, const int gSize)
   {
      computePowerQuadrature<1>(igrid, gSize);
   }
};

template <> struct EnergyHelperQuadMap<EnergySLaplR2_t>
{
   void computeQuadrature(Internal::Array& igrid, const int gSize)
   {
      computePowerQuadrature<1>(igrid, gSize);
   }
};

/// @brief Quadrature for Energy
template <> struct OpsQuadratureMap<Energy_t>
{
   using type = EnergyHelperQuadMap<Energy_t>;
};

/// @brief Quadrature for EnergyD1R1
template <> struct OpsQuadratureMap<EnergyD1R1_t>
{
   using type = EnergyHelperQuadMap<EnergyD1R1_t>;
};

/// @brief Quadrature for EnergyR2
template <> struct OpsQuadratureMap<EnergyR2_t>
{
   using type = EnergyHelperQuadMap<EnergyR2_t>;
};

/// @brief Quadrature for EnergySLaplR2
template <> struct OpsQuadratureMap<EnergySLaplR2_t>
{
   using type = EnergyHelperQuadMap<EnergySLaplR2_t>;
};

/// @brief Quadrature for Power
template <> struct OpsQuadratureMap<Power_t>
{
   using type = OpsQuadrature<Energy_t>;
};

/// @brief Quadrature for PowerD1R1
template <> struct OpsQuadratureMap<PowerD1R1_t>
{
   using type = OpsQuadrature<EnergyD1R1_t>;
};

/// @brief Quadrature for PowerR2
template <> struct OpsQuadratureMap<PowerR2_t>
{
   using type = OpsQuadrature<EnergyR2_t>;
};

/// @brief Quadrature for PowerSLaplR2
template <> struct OpsQuadratureMap<PowerSLaplR2_t>
{
   using type = OpsQuadrature<EnergySLaplR2_t>;
};

} // namespace details

} // namespace Sphere
} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC
