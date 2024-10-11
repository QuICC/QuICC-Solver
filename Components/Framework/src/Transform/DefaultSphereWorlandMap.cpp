/**
 * file DefaultSphereWorlandMap.cpp
 * @brief Source of the implementation of the Worland transform in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/DefaultSphereWorlandMap.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/P.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/SphLapl.hpp"

#include "QuICC/Transform/Poly/Worland/Integrator/P.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/P_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/R1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I4DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I4DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I2_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I2.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I2DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I2DivR1D1R1_Zero.hpp"

#include "QuICC/Transform/Poly/Worland/Reductor/EnergySLaplR2.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Energy.hpp"

#include "QuICC/Transform/Poly/Worland/Reductor/PowerSLaplR2.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/Power.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/RadialPower.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/RadialPowerDivR1.hpp"
#include "QuICC/Transform/Poly/Worland/Reductor/RadialPowerDivR1D1R1.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1R1.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/Slapl.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/I2P.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/Pol.hpp"
#include "QuICC/Transform/Forward/Q.hpp"
#include "QuICC/Transform/Forward/S.hpp"
#include "QuICC/Transform/Forward/T.hpp"
#include "QuICC/Transform/Forward/I4Q.hpp"
#include "QuICC/Transform/Forward/I4S.hpp"
#include "QuICC/Transform/Forward/I2Q.hpp"
#include "QuICC/Transform/Forward/I2S.hpp"
#include "QuICC/Transform/Forward/I2T.hpp"

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

   void DefaultSphereWorlandMap::operator()(MapType& m) const
   {
#ifdef QUICC_HAS_CUDA_BACKEND
      using backend_t = Poly::Worland::viewGpu_t;
#else
      using backend_t = Poly::Worland::viewCpu_t;
#endif

      // Create projectors
      this->addOperator<Poly::Worland::Projector::P<backend_t>>(m, Backward::P::id());
      this->addOperator<Poly::Worland::Projector::DivR1_Zero<backend_t>>(m, Backward::Overr1::id());
      this->addOperator<Poly::Worland::Projector::D1<backend_t>>(m, Backward::D1::id());
      this->addOperator<Poly::Worland::Projector::D1R1<backend_t>>(m, Backward::D1R1::id());
      this->addOperator<Poly::Worland::Projector::DivR1D1R1_Zero<backend_t>>(m, Backward::Overr1D1R1::id());
      this->addOperator<Poly::Worland::Projector::SphLapl<backend_t>>(m, Backward::Slapl::id());

      // Create integrators
      this->addOperator<Poly::Worland::Integrator::P<backend_t>>(m, Forward::P::id());
      this->addOperator<Poly::Worland::Integrator::R1_Zero<backend_t>>(m, Forward::Pol::id());
      this->addOperator<Poly::Worland::Integrator::DivR1_Zero<backend_t>>(m, Forward::Q::id());
      this->addOperator<Poly::Worland::Integrator::DivR1D1R1_Zero<backend_t>>(m, Forward::S::id());
      this->addOperator<Poly::Worland::Integrator::P_Zero<backend_t>>(m, Forward::T::id());
      this->addOperator<Poly::Worland::Integrator::I2<backend_t>>(m, Forward::I2P::id());
      this->addOperator<Poly::Worland::Integrator::I2DivR1_Zero<backend_t>>(m, Forward::I2Q::id());
      this->addOperator<Poly::Worland::Integrator::I2DivR1D1R1_Zero<backend_t>>(m, Forward::I2S::id());
      this->addOperator<Poly::Worland::Integrator::I2_Zero<backend_t>>(m, Forward::I2T::id());
      this->addOperator<Poly::Worland::Integrator::I4DivR1_Zero<backend_t>>(m, Forward::I4Q::id());
      this->addOperator<Poly::Worland::Integrator::I4DivR1D1R1_Zero<backend_t>>(m, Forward::I4S::id());

      using backendRed_t = Poly::Worland::base_t;
      this->addOperator<Poly::Worland::Reductor::EnergySLaplR2<backendRed_t>>(m, Reductor::EnergySlaplR2::id());
      this->addOperator<Poly::Worland::Reductor::EnergyD1R1<backendRed_t>>(m, Reductor::EnergyD1R1::id());
      this->addOperator<Poly::Worland::Reductor::EnergyR2<backendRed_t>>(m, Reductor::EnergyR2::id());
      this->addOperator<Poly::Worland::Reductor::Energy<backendRed_t>>(m, Reductor::Energy::id());
      this->addOperator<Poly::Worland::Reductor::PowerSLaplR2<backendRed_t>>(m, Reductor::PowerSlaplR2::id());
      this->addOperator<Poly::Worland::Reductor::PowerD1R1<backendRed_t>>(m, Reductor::PowerD1R1::id());
      this->addOperator<Poly::Worland::Reductor::PowerR2<backendRed_t>>(m, Reductor::PowerR2::id());
      this->addOperator<Poly::Worland::Reductor::Power<backendRed_t>>(m, Reductor::Power::id());
      this->addOperator<Poly::Worland::Reductor::RadialPower<backendRed_t>>(m, Reductor::RadialPower::id());
      this->addOperator<Poly::Worland::Reductor::RadialPowerDivR1<backendRed_t>>(m, Reductor::RadialPowerOverr1::id());
      this->addOperator<Poly::Worland::Reductor::RadialPowerDivR1D1R1<backendRed_t>>(m, Reductor::RadialPowerOverr1D1R1::id());
   }

}
}
