/**
 * @file OpsBuilder.hpp
 * @brief Mapping Generic Finite Differences operator builders to specific Ops
 */

#pragma once

// External includes
//

// Project includes
//
#include "SparseOp/FiniteDiff/Operator.hpp"
#include "SparseOp/FiniteDiff/OperatorWithMean.hpp"
#include "FiniteDiff/Sphere/Id.hpp"
#include "FiniteDiff/Sphere/D1.hpp"
#include "FiniteDiff/Sphere/R1.hpp"
#include "FiniteDiff/Sphere/Overr1.hpp"
#include "FiniteDiff/Sphere/Overr1D1R1.hpp"
#include "FiniteDiff/Sphere/SLapl.hpp"
#include "ViewOps/FiniteDiff/Sphere/Builder.hpp"
#include "ViewOps/FiniteDiff/Sphere/Tags.hpp"

namespace QuICC {
namespace Transform {
namespace FiniteDiff {
namespace Sphere {

/// @brief This namespace hides implementation details
namespace details {

/// @brief Generic mapping
/// @tparam VOP operator view type
/// @tparam TAG kind
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class TAG, class DIR> struct OpsBuilderMap
{
   using type = void;
};

} // namespace details

/// @brief Convenience wrapper
/// (you cannot specialize type aliases)
/// @tparam VOP operator view type
/// @tparam TAG kind
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class TAG, class DIR>
using OpsBuilder = typename details::OpsBuilderMap<VOP, TAG, DIR>::type;

// Actual mapping
namespace details {

using namespace QuICC::FiniteDiff::Sphere;

/// @brief P Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, P_t, DIR>
{
   using type =
      FiniteDiff::Sphere::Builder<VOP, QuICC::SparseOp::FiniteDiff::Operator<Id>, DIR>;
};

/// @brief D1 Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, D1_t, DIR>
{
   using type =
      FiniteDiff::Sphere::Builder<VOP, QuICC::SparseOp::FiniteDiff::Operator<D1>, DIR>;
};

/// @brief DivR1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, DivR1_t, bwd_t>
{
   using type = FiniteDiff::Sphere::Builder<VOP,
      QuICC::SparseOp::FiniteDiff::Operator<Overr1>, bwd_t>;
};

/// @brief DivR1_Zero Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, DivR1_Zero_t, DIR>
{
   using type = FiniteDiff::Sphere::Builder<VOP,
      QuICC::SparseOp::FiniteDiff::OperatorWithMean<Overr1, void>,
      DIR>;
};

/// @brief DivR1D1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, DivR1D1R1_t, bwd_t>
{
   using type = FiniteDiff::Sphere::Builder<VOP,
      QuICC::SparseOp::FiniteDiff::Operator<Overr1D1R1>, bwd_t>;
};

/// @brief DivR1D1R1_Zero Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, DivR1D1R1_Zero_t, DIR>
{
   using type = FiniteDiff::Sphere::Builder<VOP,
      QuICC::SparseOp::FiniteDiff::OperatorWithMean<Overr1D1R1, void>,
      DIR>;
};

/// @brief SphLapl Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, SphLapl_t, bwd_t>
{
   using type =
      FiniteDiff::Sphere::Builder<VOP, QuICC::SparseOp::FiniteDiff::Operator<SLapl>, bwd_t>;
};

/// @brief P_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, P_Zero_t, fwd_t>
{
   using type = FiniteDiff::Sphere::Builder<VOP,
      QuICC::SparseOp::FiniteDiff::OperatorWithMean<Id, void>, fwd_t>;
};

/// @brief R1_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, R1_Zero_t, fwd_t>
{
   using type = FiniteDiff::Sphere::Builder<VOP,
      QuICC::SparseOp::FiniteDiff::OperatorWithMean<R1, void>, fwd_t>;
};


/// @brief Generic helper in order to avoid having to pass
/// extra parameters to energy integrators
template <class VOP, class OP> struct EnergyHelperMap;

/// @brief Helper for Energy
/// It is needed in order to avoid having to pass
/// extra parameters to energy integrators
/// @tparam VOP
template <class VOP> struct EnergyHelperMap<VOP, Energy_t>
{
   void compute(VOP opView, const Internal::Array& grid)
   {
      throw std::logic_error("Not implemented");
#if 0
      Wnl fWnl(Polynomial::FiniteDiff::worland_sphenergy_t::ALPHA,
         Polynomial::FiniteDiff::worland_sphenergy_t::DBETA, -1);
      QuICC::SparseOp::FiniteDiff::OperatorWithMean<Wnl, void> denseBuilder(fWnl);
      Builder<VOP, QuICC::SparseOp::FiniteDiff::OperatorWithMean<Wnl, void>, fwd_t>
         tBuilderFwd(denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
#endif
   }
};

/// @brief Energy Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, Energy_t, fwd_t>
{
   using type = EnergyHelperMap<VOP, Energy_t>;
};

/// @brief Energy Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, Energy_t, bwd_t>
{
   using type = Builder<VOP,
      QuICC::SparseOp::FiniteDiff::Operator<
         Overr1>,
      bwd_t>;
};

/// @brief Helper for EnergyD1R1
/// It is needed in order to avoid having to pass
/// extra parameters to energy integrators
/// @tparam VOP
template <class VOP> struct EnergyHelperMap<VOP, EnergyD1R1_t>
{
   void compute(VOP opView, const Internal::Array& grid)
   {
      throw std::logic_error("Not implemented");
#if 0
      Wnl fWnl(Polynomial::FiniteDiff::worland_sphenergy_t::ALPHA,
         Polynomial::FiniteDiff::worland_sphenergy_t::DBETA, -1);
      QuICC::SparseOp::FiniteDiff::OperatorWithMean<Wnl, void> denseBuilder(fWnl);
      Builder<VOP, QuICC::SparseOp::FiniteDiff::OperatorWithMean<Wnl, void>, fwd_t>
         tBuilderFwd(denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
#endif
   }
};

/// @brief EnergyD1R1 Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, EnergyD1R1_t, fwd_t>
{
   using type = EnergyHelperMap<VOP, EnergyD1R1_t>;
};

/// @brief EnergyD1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, EnergyD1R1_t, bwd_t>
{
   using type = Builder<VOP,
      QuICC::SparseOp::FiniteDiff::Operator<
         Overr1D1R1>,
      bwd_t>;
};

/// @brief Helper for EnergyR2
/// It is needed in order to avoid having to pass
/// extra parameters to energy integrators
/// @tparam VOP
template <class VOP> struct EnergyHelperMap<VOP, EnergyR2_t>
{
   void compute(VOP opView, const Internal::Array& grid)
   {
      throw std::logic_error("Not implemented");
#if 0
      Wnl fWnl(Polynomial::FiniteDiff::worland_sphenergy_t::ALPHA,
         Polynomial::FiniteDiff::worland_sphenergy_t::DBETA, 0);
      QuICC::SparseOp::FiniteDiff::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::SparseOp::FiniteDiff::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
#endif
   }
};

/// @brief EnergyR2 Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, EnergyR2_t, fwd_t>
{
   using type = EnergyHelperMap<VOP, EnergyR2_t>;
};

/// @brief EnergyR2 Builder
/// Projector only, regular projector
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, EnergyR2_t, bwd_t>
{
   using type = OpsBuilder<VOP, P_t, bwd_t>;
};

/// @brief Helper for EnergySLaplR2
/// It is needed in order to avoid having to pass
/// extra parameters to energy integrators
/// @tparam VOP
template <class VOP> struct EnergyHelperMap<VOP, EnergySLaplR2_t>
{
   void compute(VOP opView, const Internal::Array& grid)
   {
      throw std::logic_error("Not implemented");
#if 0
      Wnl fWnl(Polynomial::FiniteDiff::worland_sphenergy_t::ALPHA,
         Polynomial::FiniteDiff::worland_sphenergy_t::DBETA, 0);
      QuICC::SparseOp::FiniteDiff::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::SparseOp::FiniteDiff::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
#endif
   }
};

/// @brief EnergySLaplR2 Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, EnergySLaplR2_t, fwd_t>
{
   using type = EnergyHelperMap<VOP, EnergySLaplR2_t>;
};

/// @brief EnergySLaplR2 Builder
/// Projector only, regular projector
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, EnergySLaplR2_t, bwd_t>
{
   using type = OpsBuilder<VOP, SphLapl_t, bwd_t>;
};

/// @brief Power Builder
/// Same setup as Energy ops
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, Power_t, DIR>
{
   using type = OpsBuilder<VOP, Energy_t, DIR>;
};

/// @brief PowerR2 Builder
/// Same setup as Energy ops
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, PowerR2_t, DIR>
{
   using type = OpsBuilder<VOP, EnergyR2_t, DIR>;
};

/// @brief PowerD1R1 Builder
/// Same setup as Energy ops
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, PowerD1R1_t, DIR>
{
   using type = OpsBuilder<VOP, EnergyD1R1_t, DIR>;
};

/// @brief PowerSLaplR2 Builder
/// Same setup as Energy ops
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, PowerSLaplR2_t, DIR>
{
   using type = OpsBuilder<VOP, EnergySLaplR2_t, DIR>;
};

/// @brief Helper for RadialPower
/// It is needed in order to avoid having to pass
/// extra parameters to energy integrators
/// @tparam VOP
template <class VOP> struct EnergyHelperMap<VOP, RadialPower_t>
{
   void compute(VOP opView, const Internal::Array& grid)
   {
      throw std::logic_error("Not implemented");
#if 0
      Wnl fWnl(Polynomial::FiniteDiff::worland_sphenergy_t::ALPHA,
         Polynomial::FiniteDiff::worland_sphenergy_t::DBETA, 0);
      QuICC::SparseOp::FiniteDiff::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::SparseOp::FiniteDiff::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
#endif
   }
};

/// @brief RadialPower Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, RadialPower_t, fwd_t>
{
   using type = EnergyHelperMap<VOP, RadialPower_t>;
};

/// @brief RadialPower Builder
/// Projector only, regular projector
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, RadialPower_t, bwd_t>
{
   using type = OpsBuilder<VOP, P_t, bwd_t>;
};

/// @brief Helper for RadialPowerDivR1
/// It is needed in order to avoid having to pass
/// extra parameters to energy integrators
/// @tparam VOP
template <class VOP> struct EnergyHelperMap<VOP, RadialPowerDivR1_t>
{
   void compute(VOP opView, const Internal::Array& grid)
   {
      throw std::logic_error("Not implemented");
#if 0
      Wnl fWnl(Polynomial::FiniteDiff::worland_sphenergy_t::ALPHA,
         Polynomial::FiniteDiff::worland_sphenergy_t::DBETA, 0);
      QuICC::SparseOp::FiniteDiff::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::SparseOp::FiniteDiff::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
#endif
   }
};

/// @brief RadialPowerDivR1 Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, RadialPowerDivR1_t, fwd_t>
{
   using type = EnergyHelperMap<VOP, RadialPowerDivR1_t>;
};

/// @brief RadialPowerDivR1 Builder
/// Projector only, regular projector
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, RadialPowerDivR1_t, bwd_t>
{
   using type = OpsBuilder<VOP, DivR1_Zero_t, bwd_t>;
};

/// @brief Helper for RadialPowerDivR1D1R1
/// It is needed in order to avoid having to pass
/// extra parameters to energy integrators
/// @tparam VOP
template <class VOP> struct EnergyHelperMap<VOP, RadialPowerDivR1D1R1_t>
{
   void compute(VOP opView, const Internal::Array& grid)
   {
      throw std::logic_error("Not implemented");
#if 0
      Wnl fWnl(Polynomial::FiniteDiff::worland_sphenergy_t::ALPHA,
         Polynomial::FiniteDiff::worland_sphenergy_t::DBETA, 0);
      QuICC::SparseOp::FiniteDiff::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::SparseOp::FiniteDiff::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
#endif
   }
};

/// @brief RadialPowerDivR1D1R1 Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, RadialPowerDivR1D1R1_t, fwd_t>
{
   using type = EnergyHelperMap<VOP, RadialPowerDivR1D1R1_t>;
};

/// @brief RadialPowerDivR1D1R1 Builder
/// Projector only, regular projector
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, RadialPowerDivR1D1R1_t, bwd_t>
{
   using type = Builder<VOP,
      QuICC::SparseOp::FiniteDiff::OperatorWithMean<
         Overr1D1R1, void>,
      bwd_t>;
};

} // namespace details

} // namespace Sphere
} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC
