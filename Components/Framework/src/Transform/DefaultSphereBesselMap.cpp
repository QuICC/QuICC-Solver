/**
 * file DefaultSphereBesselMap.cpp
 * @brief Source of the implementation of the Bessel transform in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/DefaultSphereBesselMap.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/Base/Value.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/Base/Insulating.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/Base/NoSlip.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/P.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/D1.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/SphLapl.hpp"

#include "QuICC/Transform/Poly/Bessel/Integrator/Base/Value.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/Base/Insulating.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/Base/NoSlip.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/P.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/P_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/R1_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/DivR1D1R1_Zero.hpp"

#include "QuICC/Transform/Poly/Bessel/Reductor/Base/Value.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/Insulating.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/NoSlip.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/EnergySLaplR2.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Energy.hpp"

#include "QuICC/Transform/Poly/Bessel/Reductor/PowerSLaplR2.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Power.hpp"

#include "QuICC/Transform/Backward/ValueP.hpp"
#include "QuICC/Transform/Backward/ValueOverr1.hpp"
#include "QuICC/Transform/Backward/ValueD1.hpp"
#include "QuICC/Transform/Backward/ValueOverr1D1R1.hpp"
#include "QuICC/Transform/Backward/ValueSlapl.hpp"
#include "QuICC/Transform/Backward/InsulatingP.hpp"
#include "QuICC/Transform/Backward/InsulatingOverr1.hpp"
#include "QuICC/Transform/Backward/InsulatingD1.hpp"
#include "QuICC/Transform/Backward/InsulatingOverr1D1R1.hpp"
#include "QuICC/Transform/Backward/InsulatingSlapl.hpp"
#include "QuICC/Transform/Backward/NoSlipP.hpp"
#include "QuICC/Transform/Backward/NoSlipOverr1.hpp"
#include "QuICC/Transform/Backward/NoSlipD1.hpp"
#include "QuICC/Transform/Backward/NoSlipOverr1D1R1.hpp"
#include "QuICC/Transform/Backward/NoSlipSlapl.hpp"

#include "QuICC/Transform/Forward/ValueP.hpp"
#include "QuICC/Transform/Forward/ValuePol.hpp"
#include "QuICC/Transform/Forward/ValueQ.hpp"
#include "QuICC/Transform/Forward/ValueS.hpp"
#include "QuICC/Transform/Forward/ValueT.hpp"
#include "QuICC/Transform/Forward/ValueBc1Q.hpp"
#include "QuICC/Transform/Forward/ValueBc1S.hpp"
#include "QuICC/Transform/Forward/InsulatingP.hpp"
#include "QuICC/Transform/Forward/InsulatingPol.hpp"
#include "QuICC/Transform/Forward/InsulatingQ.hpp"
#include "QuICC/Transform/Forward/InsulatingS.hpp"
#include "QuICC/Transform/Forward/InsulatingT.hpp"
#include "QuICC/Transform/Forward/NoSlipP.hpp"
#include "QuICC/Transform/Forward/NoSlipPol.hpp"
#include "QuICC/Transform/Forward/NoSlipQ.hpp"
#include "QuICC/Transform/Forward/NoSlipS.hpp"
#include "QuICC/Transform/Forward/NoSlipT.hpp"
#include "QuICC/Transform/Forward/NoSlipBc1Q.hpp"
#include "QuICC/Transform/Forward/NoSlipBc1S.hpp"

#include "QuICC/Transform/Reductor/ValueEnergy.hpp"
#include "QuICC/Transform/Reductor/ValueEnergyR2.hpp"
#include "QuICC/Transform/Reductor/ValueEnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/ValueEnergySlaplR2.hpp"
#include "QuICC/Transform/Reductor/InsulatingEnergy.hpp"
#include "QuICC/Transform/Reductor/InsulatingEnergyR2.hpp"
#include "QuICC/Transform/Reductor/InsulatingEnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/InsulatingEnergySlaplR2.hpp"
#include "QuICC/Transform/Reductor/NoSlipEnergy.hpp"
#include "QuICC/Transform/Reductor/NoSlipEnergyR2.hpp"
#include "QuICC/Transform/Reductor/NoSlipEnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/NoSlipEnergySlaplR2.hpp"

#include "QuICC/Transform/Reductor/ValuePower.hpp"
#include "QuICC/Transform/Reductor/ValuePowerR2.hpp"
#include "QuICC/Transform/Reductor/ValuePowerD1R1.hpp"
#include "QuICC/Transform/Reductor/ValuePowerSlaplR2.hpp"
#include "QuICC/Transform/Reductor/InsulatingPower.hpp"
#include "QuICC/Transform/Reductor/InsulatingPowerR2.hpp"
#include "QuICC/Transform/Reductor/InsulatingPowerD1R1.hpp"
#include "QuICC/Transform/Reductor/InsulatingPowerSlaplR2.hpp"
#include "QuICC/Transform/Reductor/NoSlipPower.hpp"
#include "QuICC/Transform/Reductor/NoSlipPowerR2.hpp"
#include "QuICC/Transform/Reductor/NoSlipPowerD1R1.hpp"
#include "QuICC/Transform/Reductor/NoSlipPowerSlaplR2.hpp"

namespace QuICC {

namespace Transform {

   void DefaultSphereBesselMap::operator()(MapType& m) const
   {
#ifdef QUICC_HAS_CUDA_BACKEND
#error "GPU Backend for Bessel transform not implemented"
#else
      using backend_t = Poly::Bessel::base_t;
#endif

      // Create projectors
      namespace proj_ns = Poly::Bessel::Projector; 
      // ... Value BC
      this->addOperator<proj_ns::Value<proj_ns::P<backend_t>>>(m, Backward::ValueP::id());
      this->addOperator<proj_ns::Value<proj_ns::DivR1_Zero<backend_t>>>(m, Backward::ValueOverr1::id());
      this->addOperator<proj_ns::Value<proj_ns::D1<backend_t>>>(m, Backward::ValueD1::id());
      this->addOperator<proj_ns::Value<proj_ns::DivR1D1R1_Zero<backend_t>>>(m, Backward::ValueOverr1D1R1::id());
      this->addOperator<proj_ns::Value<proj_ns::SphLapl<backend_t>>>(m, Backward::ValueSlapl::id());
      // ... Insulating BC
      this->addOperator<proj_ns::Insulating<proj_ns::P<backend_t>>>(m, Backward::InsulatingP::id());
      this->addOperator<proj_ns::Insulating<proj_ns::DivR1_Zero<backend_t>>>(m, Backward::InsulatingOverr1::id());
      this->addOperator<proj_ns::Insulating<proj_ns::D1<backend_t>>>(m, Backward::InsulatingD1::id());
      this->addOperator<proj_ns::Insulating<proj_ns::DivR1D1R1_Zero<backend_t>>>(m, Backward::InsulatingOverr1D1R1::id());
      this->addOperator<proj_ns::Insulating<proj_ns::SphLapl<backend_t>>>(m, Backward::InsulatingSlapl::id());
      // ... NoSlip BC
      this->addOperator<proj_ns::NoSlip<proj_ns::P<backend_t>>>(m, Backward::NoSlipP::id());
      this->addOperator<proj_ns::NoSlip<proj_ns::DivR1_Zero<backend_t>>>(m, Backward::NoSlipOverr1::id());
      this->addOperator<proj_ns::NoSlip<proj_ns::D1<backend_t>>>(m, Backward::NoSlipD1::id());
      this->addOperator<proj_ns::NoSlip<proj_ns::DivR1D1R1_Zero<backend_t>>>(m, Backward::NoSlipOverr1D1R1::id());
      this->addOperator<proj_ns::NoSlip<proj_ns::SphLapl<backend_t>>>(m, Backward::NoSlipSlapl::id());

      // Create integrators
      namespace intg_ns = Poly::Bessel::Integrator; 
      // ... Value BC
      this->addOperator<intg_ns::Value<intg_ns::P<backend_t>>>(m, Forward::ValueP::id());
      this->addOperator<intg_ns::Value<intg_ns::R1_Zero<backend_t>>>(m, Forward::ValuePol::id());
      this->addOperator<intg_ns::Value<intg_ns::DivR1_Zero<backend_t>>>(m, Forward::ValueQ::id());
      this->addOperator<intg_ns::Value<intg_ns::DivR1_Zero<backend_t>,1>>(m, Forward::ValueBc1Q::id());
      this->addOperator<intg_ns::Value<intg_ns::DivR1D1R1_Zero<backend_t>>>(m, Forward::ValueS::id());
      this->addOperator<intg_ns::Value<intg_ns::DivR1D1R1_Zero<backend_t>,1>>(m, Forward::ValueBc1S::id());
      this->addOperator<intg_ns::Value<intg_ns::P_Zero<backend_t>>>(m, Forward::ValueT::id());
      // ... Insulating BC
      this->addOperator<intg_ns::Insulating<intg_ns::P<backend_t>>>(m, Forward::InsulatingP::id());
      this->addOperator<intg_ns::Insulating<intg_ns::R1_Zero<backend_t>>>(m, Forward::InsulatingPol::id());
      this->addOperator<intg_ns::Insulating<intg_ns::DivR1_Zero<backend_t>>>(m, Forward::InsulatingQ::id());
      this->addOperator<intg_ns::Insulating<intg_ns::DivR1D1R1_Zero<backend_t>>>(m, Forward::InsulatingS::id());
      this->addOperator<intg_ns::Insulating<intg_ns::P_Zero<backend_t>>>(m, Forward::InsulatingT::id());
      // ... NoSlip BC
      this->addOperator<intg_ns::NoSlip<intg_ns::P<backend_t>>>(m, Forward::NoSlipP::id());
      this->addOperator<intg_ns::NoSlip<intg_ns::R1_Zero<backend_t>>>(m, Forward::NoSlipPol::id());
      this->addOperator<intg_ns::NoSlip<intg_ns::DivR1_Zero<backend_t>>>(m, Forward::NoSlipQ::id());
      this->addOperator<intg_ns::NoSlip<intg_ns::DivR1_Zero<backend_t>,0,1>>(m, Forward::NoSlipBc1Q::id());
      this->addOperator<intg_ns::NoSlip<intg_ns::DivR1D1R1_Zero<backend_t>>>(m, Forward::NoSlipS::id());
      this->addOperator<intg_ns::NoSlip<intg_ns::DivR1D1R1_Zero<backend_t>,0,1>>(m, Forward::NoSlipBc1S::id());
      this->addOperator<intg_ns::NoSlip<intg_ns::P_Zero<backend_t>>>(m, Forward::NoSlipT::id());

      // Create reductors
      namespace red_ns = Poly::Bessel::Reductor; 
      // ... Energy
      // ... ... Value BC
      this->addOperator<red_ns::Value<red_ns::EnergySLaplR2<backend_t>>>(m, Reductor::ValueEnergySlaplR2::id());
      this->addOperator<red_ns::Value<red_ns::EnergyD1R1<backend_t>>>(m, Reductor::ValueEnergyD1R1::id());
      this->addOperator<red_ns::Value<red_ns::EnergyR2<backend_t>>>(m, Reductor::ValueEnergyR2::id());
      this->addOperator<red_ns::Value<red_ns::Energy<backend_t>>>(m, Reductor::ValueEnergy::id());
      // ... ... Insulating BC
      this->addOperator<red_ns::Insulating<red_ns::EnergySLaplR2<backend_t>>>(m, Reductor::InsulatingEnergySlaplR2::id());
      this->addOperator<red_ns::Insulating<red_ns::EnergyD1R1<backend_t>>>(m, Reductor::InsulatingEnergyD1R1::id());
      this->addOperator<red_ns::Insulating<red_ns::EnergyR2<backend_t>>>(m, Reductor::InsulatingEnergyR2::id());
      this->addOperator<red_ns::Insulating<red_ns::Energy<backend_t>>>(m, Reductor::InsulatingEnergy::id());
      // ... ... NoSlip BC
      this->addOperator<red_ns::NoSlip<red_ns::EnergySLaplR2<backend_t>>>(m, Reductor::NoSlipEnergySlaplR2::id());
      this->addOperator<red_ns::NoSlip<red_ns::EnergyD1R1<backend_t>>>(m, Reductor::NoSlipEnergyD1R1::id());
      this->addOperator<red_ns::NoSlip<red_ns::EnergyR2<backend_t>>>(m, Reductor::NoSlipEnergyR2::id());
      this->addOperator<red_ns::NoSlip<red_ns::Energy<backend_t>>>(m, Reductor::NoSlipEnergy::id());
      // ... Power
      // ... ... Value BC
      this->addOperator<red_ns::Value<red_ns::PowerSLaplR2<backend_t>>>(m, Reductor::ValuePowerSlaplR2::id());
      this->addOperator<red_ns::Value<red_ns::PowerD1R1<backend_t>>>(m, Reductor::ValuePowerD1R1::id());
      this->addOperator<red_ns::Value<red_ns::PowerR2<backend_t>>>(m, Reductor::ValuePowerR2::id());
      this->addOperator<red_ns::Value<red_ns::Power<backend_t>>>(m, Reductor::ValuePower::id());
      // ... ... Insulating BC
      this->addOperator<red_ns::Insulating<red_ns::PowerSLaplR2<backend_t>>>(m, Reductor::InsulatingPowerSlaplR2::id());
      this->addOperator<red_ns::Insulating<red_ns::PowerD1R1<backend_t>>>(m, Reductor::InsulatingPowerD1R1::id());
      this->addOperator<red_ns::Insulating<red_ns::PowerR2<backend_t>>>(m, Reductor::InsulatingPowerR2::id());
      this->addOperator<red_ns::Insulating<red_ns::Power<backend_t>>>(m, Reductor::InsulatingPower::id());
      // ... ... NoSlip BC
      this->addOperator<red_ns::NoSlip<red_ns::PowerSLaplR2<backend_t>>>(m, Reductor::NoSlipPowerSlaplR2::id());
      this->addOperator<red_ns::NoSlip<red_ns::PowerD1R1<backend_t>>>(m, Reductor::NoSlipPowerD1R1::id());
      this->addOperator<red_ns::NoSlip<red_ns::PowerR2<backend_t>>>(m, Reductor::NoSlipPowerR2::id());
      this->addOperator<red_ns::NoSlip<red_ns::Power<backend_t>>>(m, Reductor::NoSlipPower::id());
   }

}
}
