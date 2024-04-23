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
#include "QuICC/Transform/Poly/Bessel/Projector/P.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/D1.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Projector/SphLapl.hpp"

#include "QuICC/Transform/Poly/Bessel/Integrator/Base/Value.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/Base/Insulating.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/P.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/P_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/R1_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Bessel/Integrator/DivR1D1R1_Zero.hpp"

#include "QuICC/Transform/Poly/Bessel/Reductor/Base/Value.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Base/Insulating.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/EnergySLaplR2.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Energy.hpp"

#include "QuICC/Transform/Poly/Bessel/Reductor/PowerSLaplR2.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Poly/Bessel/Reductor/Power.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/Slapl.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/Pol.hpp"
#include "QuICC/Transform/Forward/Q.hpp"
#include "QuICC/Transform/Forward/S.hpp"
#include "QuICC/Transform/Forward/T.hpp"

#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/EnergySlaplR2.hpp"

#include "QuICC/Transform/Reductor/Power.hpp"
#include "QuICC/Transform/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/Reductor/PowerSlaplR2.hpp"

#include "QuICC/Transform/Reductor/RadialPower.hpp"
#include "QuICC/Transform/Reductor/RadialPowerOverr1.hpp"
#include "QuICC/Transform/Reductor/RadialPowerOverr1D1R1.hpp"

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
      this->addOperator<proj_ns::Value<proj_ns::P<backend_t>>>(m, Backward::P::id());
      this->addOperator<proj_ns::Value<proj_ns::DivR1_Zero<backend_t>>>(m, Backward::Overr1::id());
      this->addOperator<proj_ns::Value<proj_ns::D1<backend_t>>>(m, Backward::D1::id());
      this->addOperator<proj_ns::Value<proj_ns::DivR1D1R1_Zero<backend_t>>>(m, Backward::Overr1D1R1::id());
      this->addOperator<proj_ns::Value<proj_ns::SphLapl<backend_t>>>(m, Backward::Slapl::id());

      // Create integrators
      namespace intg_ns = Poly::Bessel::Integrator; 
      this->addOperator<intg_ns::Value<intg_ns::P<backend_t>>>(m, Forward::P::id());
      this->addOperator<intg_ns::Value<intg_ns::R1_Zero<backend_t>>>(m, Forward::Pol::id());
      this->addOperator<intg_ns::Value<intg_ns::DivR1_Zero<backend_t>>>(m, Forward::Q::id());
      this->addOperator<intg_ns::Value<intg_ns::DivR1D1R1_Zero<backend_t>>>(m, Forward::S::id());
      this->addOperator<intg_ns::Value<intg_ns::P_Zero<backend_t>>>(m, Forward::T::id());

      // Create reductors
      namespace red_ns = Poly::Bessel::Reductor; 
      this->addOperator<red_ns::Value<red_ns::EnergySLaplR2<backend_t>>>(m, Reductor::EnergySlaplR2::id());
      this->addOperator<red_ns::Value<red_ns::EnergyD1R1<backend_t>>>(m, Reductor::EnergyD1R1::id());
      this->addOperator<red_ns::Value<red_ns::EnergyR2<backend_t>>>(m, Reductor::EnergyR2::id());
      this->addOperator<red_ns::Value<red_ns::Energy<backend_t>>>(m, Reductor::Energy::id());
      this->addOperator<red_ns::Value<red_ns::PowerSLaplR2<backend_t>>>(m, Reductor::PowerSlaplR2::id());
      this->addOperator<red_ns::Value<red_ns::PowerD1R1<backend_t>>>(m, Reductor::PowerD1R1::id());
      this->addOperator<red_ns::Value<red_ns::PowerR2<backend_t>>>(m, Reductor::PowerR2::id());
      this->addOperator<red_ns::Value<red_ns::Power<backend_t>>>(m, Reductor::Power::id());
   }

}
}
