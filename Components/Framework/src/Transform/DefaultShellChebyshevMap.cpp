/**
 * file DefaultShellChebyshevMap.cpp
 * @brief Source of the implementation of the Chebyshev transform in a spherical
 * shell
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1R1.hpp"
#include "QuICC/Transform/Backward/D2.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/Overr1D0.hpp" // equivalent to the above
#include "QuICC/Transform/Backward/Overr1D1.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/Overr2D1R1.hpp" 
#include "QuICC/Transform/Backward/OverrSqD1R1.hpp" // equivalent to the above
#include "QuICC/Transform/Backward/Overr2.hpp"
#include "QuICC/Transform/Backward/OverrSq.hpp" // equivalent to the above
#include "QuICC/Transform/Backward/Overr1D2R1.hpp"
#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Slaplr.hpp"
#include "QuICC/Transform/DefaultShellChebyshevMap.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y1D1Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y2D1Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y2_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4Y3D1Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4Y3_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1D1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY2D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1D2Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/SphRadLapl.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Energy.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyD1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergySLaplR2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyY2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/RadialPower.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/RadialPowerDivY1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/RadialPowerDivY1D1Y1.hpp"
#include "QuICC/Transform/Forward/I2Q.hpp"
#include "QuICC/Transform/Forward/I2S.hpp"
#include "QuICC/Transform/Forward/I2T.hpp"
#include "QuICC/Transform/Forward/I2rQ.hpp"
#include "QuICC/Transform/Forward/I2rS.hpp"
#include "QuICC/Transform/Forward/I4Q.hpp"
#include "QuICC/Transform/Forward/I4S.hpp"
#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/Pol.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/T.hpp"
#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"
#include "QuICC/Transform/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Reductor/EnergySlaplR2.hpp"
#include "QuICC/Transform/Reductor/RadialPower.hpp"
#include "QuICC/Transform/Reductor/RadialPowerOverr1.hpp"
#include "QuICC/Transform/Reductor/RadialPowerOverr1D1R1.hpp"

#include "QuICC/Transform/Reductor/PowerR2.hpp"

#include "QuICC/Transform/Reductor/Spectrum.hpp"

namespace QuICC {

namespace Transform {

void DefaultShellChebyshevMap::operator()(MapType& m) const
{
#ifdef QUICC_HAS_CUDA_BACKEND
   // using backend_t = Fft::Chebyshev::LinearMap::viewGpu_t;
   using backend_t = Fft::Chebyshev::LinearMap::base_t;
#else
   // using backend_t = Fft::Chebyshev::LinearMap::viewCpu_t;
   using backend_t = Fft::Chebyshev::LinearMap::base_t;
#endif

   // Create projectors
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::P<backend_t>>(m,
      Backward::P::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1<backend_t>>(m,
      Backward::Overr1::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1<backend_t>>(m,
      Backward::Overr1D0::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY2<backend_t>>(m,
      Backward::Overr2::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY2<backend_t>>(m,
      Backward::OverrSq::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::D1<backend_t>>(m,
      Backward::D1::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::D1Y1<backend_t>>(m,
      Backward::D1R1::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::D2<backend_t>>(m,
      Backward::D2::id());
   this
      ->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1D1Y1<backend_t>>(
         m, Backward::Overr1D1R1::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1D1<backend_t>>(m, 
      Backward::Overr1D1::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY2D1Y1<backend_t>>(m, 
      Backward::Overr2D1R1::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Projector::SphRadLapl<backend_t>>(m,
      Backward::Slaplr::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1D2Y1<backend_t>>(m, 
      Backward::Overr1D2R1::id());

   // Create integrators
   this->addOperator<Fft::Chebyshev::LinearMap::Integrator::P<backend_t>>(m,
      Forward::P::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Integrator::Y1_Zero<backend_t>>(
      m, Forward::Pol::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Integrator::P_Zero<backend_t>>(
      m, Forward::T::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Integrator::I4Y3_Zero<backend_t>>(m,
      Forward::I4Q::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Integrator::I4Y3D1Y1_Zero<backend_t>>(m,
      Forward::I4S::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Integrator::I2Y2_Zero<backend_t>>(m,
      Forward::I2T::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Integrator::I2Y1_Zero<backend_t>>(m,
      Forward::I2Q::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Integrator::I2Y1D1Y1_Zero<backend_t>>(m,
      Forward::I2S::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Integrator::I2Y2_Zero<backend_t>>(m,
      Forward::I2rQ::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Integrator::I2Y2D1Y1_Zero<backend_t>>(m,
      Forward::I2rS::id());

   // Create reductors
   this->addOperator<Fft::Chebyshev::LinearMap::Reductor::Energy<backend_t>>(m,
      Reductor::Energy::id());
   this
      ->addOperator<Fft::Chebyshev::LinearMap::Reductor::EnergyD1Y1<backend_t>>(
         m, Reductor::EnergyD1R1::id());
   this->addOperator<Fft::Chebyshev::LinearMap::Reductor::EnergyY2<backend_t>>(
      m, Reductor::EnergyR2::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Reductor::EnergySLaplR2<backend_t>>(m,
      Reductor::EnergySlaplR2::id());

   this->addOperator<
      Fft::Chebyshev::LinearMap::Reductor::RadialPower<backend_t>>(m,
      Reductor::RadialPower::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Reductor::RadialPowerDivY1<backend_t>>(m,
      Reductor::RadialPowerOverr1::id());
   this->addOperator<
      Fft::Chebyshev::LinearMap::Reductor::RadialPowerDivY1D1Y1<backend_t>>(m,
      Reductor::RadialPowerOverr1D1R1::id());
}

} // namespace Transform
} // namespace QuICC
