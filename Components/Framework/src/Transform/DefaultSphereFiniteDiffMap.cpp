/**
 * file DefaultSphereFiniteDiffMap.cpp
 * @brief Source of the implementation of the Finite Differences transform in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1R1.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Slapl.hpp"
#include "QuICC/Transform/DefaultSphereFiniteDiffMap.hpp"
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/Pol.hpp"
#include "QuICC/Transform/Forward/Q.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/S.hpp"
#include "QuICC/Transform/Forward/T.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Integrator/DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Integrator/DivR1_Zero.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Integrator/P.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Integrator/P_Zero.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Integrator/R1_Zero.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Projector/D1.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Projector/DivR1D1R1_Zero.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Projector/DivR1_Zero.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Projector/P.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Projector/SphLapl.hpp"
#ifdef FD_ENERGY_IMPLEMENTED
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/Energy.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/EnergySLaplR2.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/Power.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/PowerR2.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/PowerSLaplR2.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/RadialPower.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/RadialPowerDivR1.hpp"
#include "QuICC/Transform/FiniteDiff/Sphere/Reductor/RadialPowerDivR1D1R1.hpp"
#endif
#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Reductor/EnergySlaplR2.hpp"
#include "QuICC/Transform/Reductor/Power.hpp"
#include "QuICC/Transform/Reductor/PowerD1R1.hpp"
#include "QuICC/Transform/Reductor/PowerR2.hpp"
#include "QuICC/Transform/Reductor/PowerSlaplR2.hpp"
#include "QuICC/Transform/Reductor/RadialPower.hpp"
#include "QuICC/Transform/Reductor/RadialPowerOverr1.hpp"
#include "QuICC/Transform/Reductor/RadialPowerOverr1D1R1.hpp"

namespace QuICC {

namespace Transform {

void DefaultSphereFiniteDiffMap::operator()(MapType& m) const
{
#ifdef QUICC_HAS_CUDA_BACKEND
   using backend_t = FiniteDiff::Sphere::viewGpu_t;
#else
   using backend_t = FiniteDiff::Sphere::viewCpu_t;
#endif

   // Create projectors
   this->addOperator<FiniteDiff::Sphere::Projector::P<backend_t>>(m,
      Backward::P::id());
   this->addOperator<FiniteDiff::Sphere::Projector::DivR1_Zero<backend_t>>(m,
      Backward::Overr1::id());
   this->addOperator<FiniteDiff::Sphere::Projector::D1<backend_t>>(m,
      Backward::D1::id());
   this->addOperator<FiniteDiff::Sphere::Projector::DivR1D1R1_Zero<backend_t>>(m,
      Backward::Overr1D1R1::id());
   this->addOperator<FiniteDiff::Sphere::Projector::SphLapl<backend_t>>(m,
      Backward::Slapl::id());

   // Create integrators
   this->addOperator<FiniteDiff::Sphere::Integrator::P<backend_t>>(m,
      Forward::P::id());
   this->addOperator<FiniteDiff::Sphere::Integrator::R1_Zero<backend_t>>(m,
      Forward::Pol::id());
   this->addOperator<FiniteDiff::Sphere::Integrator::DivR1_Zero<backend_t>>(m,
      Forward::Q::id());
   this->addOperator<FiniteDiff::Sphere::Integrator::DivR1D1R1_Zero<backend_t>>(m,
      Forward::S::id());
   this->addOperator<FiniteDiff::Sphere::Integrator::P_Zero<backend_t>>(m,
      Forward::T::id());

   // Create reductors
#ifdef FD_ENERGY_IMPLEMENTED
   this->addOperator<FiniteDiff::Sphere::Reductor::EnergySLaplR2<backend_t>>(m,
      Reductor::EnergySlaplR2::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::EnergyD1R1<backend_t>>(m,
      Reductor::EnergyD1R1::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::EnergyR2<backend_t>>(m,
      Reductor::EnergyR2::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::Energy<backend_t>>(m,
      Reductor::Energy::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::PowerSLaplR2<backend_t>>(m,
      Reductor::PowerSlaplR2::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::PowerD1R1<backend_t>>(m,
      Reductor::PowerD1R1::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::PowerR2<backend_t>>(m,
      Reductor::PowerR2::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::Power<backend_t>>(m,
      Reductor::Power::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::RadialPower<backend_t>>(m,
      Reductor::RadialPower::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::RadialPowerDivR1<backend_t>>(m,
      Reductor::RadialPowerOverr1::id());
   this->addOperator<FiniteDiff::Sphere::Reductor::RadialPowerDivR1D1R1<backend_t>>(
      m, Reductor::RadialPowerOverr1D1R1::id());
#endif
}

} // namespace Transform
} // namespace QuICC
