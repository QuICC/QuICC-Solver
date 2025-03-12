/**
 * file DefaultShellChebyshevMap.cpp
 * @brief Source of the implementation of the Chebyshev transform in a spherical shell
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/DefaultShellChebyshevMap.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/SphRadLapl.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
                                              
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4Y3_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4Y3D1Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y2_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y1D1Y1_Zero.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2Y2D1Y1_Zero.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Energy.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyY2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyD1Y1.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/R1.hpp"
#include "QuICC/Transform/Forward/Pol.hpp"
#include "QuICC/Transform/Forward/I4Q.hpp"
#include "QuICC/Transform/Forward/I4S.hpp"
#include "QuICC/Transform/Forward/I2Q.hpp"
#include "QuICC/Transform/Forward/I2S.hpp"
#include "QuICC/Transform/Forward/I2rQ.hpp"
#include "QuICC/Transform/Forward/I2rS.hpp"
#include "QuICC/Transform/Forward/I2T.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/Overr2.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1R1.hpp"
#include "QuICC/Transform/Backward/D2.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/Slaplr.hpp"

#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyR2.hpp"
#include "QuICC/Transform/Reductor/EnergyD1R1.hpp"

namespace QuICC {

namespace Transform {

   void DefaultShellChebyshevMap::operator()(MapType& m) const
   {
      // Create projectors
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::P>(m, Backward::P::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1>(m, Backward::Overr1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY2>(m, Backward::Overr2::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::D1>(m, Backward::D1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::D1Y1>(m, Backward::D1R1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::D2>(m, Backward::D2::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1D1Y1>(m, Backward::Overr1D1R1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::SphRadLapl>(m, Backward::Slaplr::id());

      // Create integrators
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::P>(m, Forward::P::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::Y1>(m, Forward::Pol::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I4Y3_Zero>(m, Forward::I4Q::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I4Y3D1Y1_Zero>(m, Forward::I4S::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I2Y2_Zero>(m, Forward::I2T::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I2Y1_Zero>(m, Forward::I2Q::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I2Y1D1Y1_Zero>(m, Forward::I2S::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I2Y2_Zero>(m, Forward::I2rQ::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I2Y2D1Y1_Zero>(m, Forward::I2rS::id());
      

      // Create reductors
      this->addOperator<Fft::Chebyshev::LinearMap::Reductor::Energy>(m, Reductor::Energy::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Reductor::EnergyD1Y1>(m, Reductor::EnergyD1R1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Reductor::EnergyY2>(m, Reductor::EnergyR2::id());
   }
 
} // Transform
} // QuICC
