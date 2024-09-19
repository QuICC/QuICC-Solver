/**
 * @file OpsBuilder.hpp
 * @brief Mapping Generic Worland operator builders to specific Ops
 */

#pragma once

// External includes
//

// Project includes
//
#include "DenseSM/Worland/Operator.hpp"
#include "DenseSM/Worland/OperatorIN.hpp"
#include "DenseSM/Worland/OperatorINWithMean.hpp"
#include "DenseSM/Worland/OperatorWithMean.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Polynomial/Worland/claplhWnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/dclaplhWnl.hpp"
#include "QuICC/Polynomial/Worland/drWnl.hpp"
#include "QuICC/Polynomial/Worland/dr_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/rWnl.hpp"
#include "QuICC/Polynomial/Worland/r_1Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1claplhWnl.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "ViewOps/Worland/Builder.hpp"
#include "ViewOps/Worland/Tags.hpp"

namespace QuICC {
namespace Transform {
namespace Worland {

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

using namespace QuICC::Polynomial::Worland;

/// @brief P Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, P_t, DIR>
{
   using type =
      Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<Wnl>, DIR>;
};

/// @brief D1 Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR> struct OpsBuilderMap<VOP, D1_t, DIR>
{
   using type =
      Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<dWnl>, DIR>;
};

/// @brief D1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, D1R1_t, bwd_t>
{
   using type =
      Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<drWnl>, bwd_t>;
};

/// @brief D1_P Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, D1_P_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<dWnl, Wnl>, bwd_t>;
};

/// @brief DivR1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, DivR1_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::Operator<r_1Wnl<recurrence_t>>, bwd_t>;
};

/// @brief DivR1_Zero Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, DivR1_Zero_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<r_1Wnl<recurrence_t>, void>,
      bwd_t>;
};

/// @brief DivR1_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, DivR1_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<r_1Wnl<implicit_t>, void>,
      fwd_t>;
};

/// @brief DivR1D1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, DivR1D1R1_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::Operator<r_1drWnl<recurrence_t>>, bwd_t>;
};

/// @brief DivR1D1R1_Zero Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, DivR1D1R1_Zero_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<r_1drWnl<recurrence_t>, void>,
      bwd_t>;
};

/// @brief DivR1D1R1_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, DivR1D1R1_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<r_1drWnl<implicit_t>, void>,
      fwd_t>;
};

/// @brief SphLapl Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, SphLapl_t, bwd_t>
{
   using type =
      Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<slaplWnl>, bwd_t>;
};

/// @brief CylLaplh Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, CylLaplh_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::Operator<claplhWnl>, bwd_t>;
};

/// @brief CylLaplh_DivR1D1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, CylLaplh_DivR1D1R1_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<claplhWnl,
         r_1drWnl<recurrence_t>>,
      bwd_t>;
};

/// @brief D1CylLaplh Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, D1CylLaplh_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::Operator<dclaplhWnl>, bwd_t>;
};

/// @brief D1CylLaplh_D1DivR1D1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, D1CylLaplh_D1DivR1D1R1_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<dclaplhWnl, dr_1drWnl>, bwd_t>;
};

/// @brief DivR1CylLaplh_Zero Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, DivR1CylLaplh_Zero_t, bwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<r_1claplhWnl, void>, bwd_t>;
};

/// @brief P_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, P_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<Wnl, void>, fwd_t>;
};

/// @brief R1_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, R1_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorWithMean<rWnl, void>, fwd_t>;
};

/// @brief I2 Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, I2_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorIN<Wnl, ::QuICC::SparseSM::Worland::I2>,
      fwd_t>;
};
/// @brief I2_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, I2_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorINWithMean<Wnl,
         ::QuICC::SparseSM::Worland::I2>,
      fwd_t>;
};

/// @brief I2DivR1_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, I2DivR1_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorINWithMean<r_1Wnl<implicit_t>,
         ::QuICC::SparseSM::Worland::I2>,
      fwd_t>;
};

/// @brief I4DivR1_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, I4DivR1_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorINWithMean<r_1Wnl<implicit_t>,
         ::QuICC::SparseSM::Worland::I4>,
      fwd_t>;
};

/// @brief I6DivR1_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, I6DivR1_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorINWithMean<r_1Wnl<implicit_t>,
         ::QuICC::SparseSM::Worland::I6>,
      fwd_t>;
};

/// @brief I2DivR1D1R1_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, I2DivR1D1R1_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorINWithMean<r_1drWnl<implicit_t>,
         ::QuICC::SparseSM::Worland::I2>,
      fwd_t>;
};

/// @brief I4DivR1D1R1_Zero Builder
/// Integrator only
/// @tparam VOP operator view type
template <class VOP> struct OpsBuilderMap<VOP, I4DivR1D1R1_Zero_t, fwd_t>
{
   using type = Worland::Builder<VOP,
      QuICC::DenseSM::Worland::OperatorINWithMean<r_1drWnl<implicit_t>,
         ::QuICC::SparseSM::Worland::I4>,
      fwd_t>;
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
   void compute(VOP opView, const Internal::Array& grid,
      const Internal::Array& weights)
   {
      Wnl fWnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
         Polynomial::Worland::worland_sphenergy_t::DBETA, -1);
      QuICC::DenseSM::Worland::OperatorWithMean<Wnl, void> denseBuilder(fWnl);
      Builder<VOP, QuICC::DenseSM::Worland::OperatorWithMean<Wnl, void>, fwd_t>
         tBuilderFwd(denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
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
      QuICC::DenseSM::Worland::Operator<
         r_1Wnl<QuICC::Polynomial::Worland::recurrence_t>>,
      bwd_t>;
};

/// @brief Helper for EnergyD1R1
/// It is needed in order to avoid having to pass
/// extra parameters to energy integrators
/// @tparam VOP
template <class VOP> struct EnergyHelperMap<VOP, EnergyD1R1_t>
{
   void compute(VOP opView, const Internal::Array& grid,
      const Internal::Array& weights)
   {
      Wnl fWnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
         Polynomial::Worland::worland_sphenergy_t::DBETA, -1);
      QuICC::DenseSM::Worland::OperatorWithMean<Wnl, void> denseBuilder(fWnl);
      Builder<VOP, QuICC::DenseSM::Worland::OperatorWithMean<Wnl, void>, fwd_t>
         tBuilderFwd(denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
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
      QuICC::DenseSM::Worland::Operator<
         r_1drWnl<QuICC::Polynomial::Worland::recurrence_t>>,
      bwd_t>;
};

/// @brief Helper for EnergyR2
/// It is needed in order to avoid having to pass
/// extra parameters to energy integrators
/// @tparam VOP
template <class VOP> struct EnergyHelperMap<VOP, EnergyR2_t>
{
   void compute(VOP opView, const Internal::Array& grid,
      const Internal::Array& weights)
   {
      Wnl fWnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
         Polynomial::Worland::worland_sphenergy_t::DBETA, 0);
      QuICC::DenseSM::Worland::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
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
   void compute(VOP opView, const Internal::Array& grid,
      const Internal::Array& weights)
   {
      Wnl fWnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
         Polynomial::Worland::worland_sphenergy_t::DBETA, 0);
      QuICC::DenseSM::Worland::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
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
   void compute(VOP opView, const Internal::Array& grid,
      const Internal::Array& weights)
   {
      Wnl fWnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
         Polynomial::Worland::worland_sphenergy_t::DBETA, 0);
      QuICC::DenseSM::Worland::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
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
   void compute(VOP opView, const Internal::Array& grid,
      const Internal::Array& weights)
   {
      Wnl fWnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
         Polynomial::Worland::worland_sphenergy_t::DBETA, 0);
      QuICC::DenseSM::Worland::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
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
   void compute(VOP opView, const Internal::Array& grid,
      const Internal::Array& weights)
   {
      Wnl fWnl(Polynomial::Worland::worland_sphenergy_t::ALPHA,
         Polynomial::Worland::worland_sphenergy_t::DBETA, 0);
      QuICC::DenseSM::Worland::Operator<Wnl> denseBuilder(fWnl);
      Builder<VOP, QuICC::DenseSM::Worland::Operator<Wnl>, fwd_t> tBuilderFwd(
         denseBuilder);
      tBuilderFwd.compute(opView, grid, weights);
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
      QuICC::DenseSM::Worland::OperatorWithMean<
         r_1drWnl<QuICC::Polynomial::Worland::recurrence_t>, void>,
      bwd_t>;
};

} // namespace details

} // namespace Worland
} // namespace Transform
} // namespace QuICC
