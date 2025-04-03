/**
 * file DefaultCartesianChebyshevMap.cpp
 * @brief Source of the implementation of the Chebyshev transform in a cartesian box
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/DefaultCartesianChebyshevMap.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2D1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4D1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I2_I2D1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/I4D1_I2.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Energy.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/EnergyD1.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D2.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/I2P.hpp"
#include "QuICC/Transform/Forward/I2D1.hpp"
#include "QuICC/Transform/Forward/I2ZI2D1.hpp"
#include "QuICC/Transform/Forward/I4P.hpp"
#include "QuICC/Transform/Forward/I4D1.hpp"
#include "QuICC/Transform/Forward/I4D1ZI2.hpp"

#include "QuICC/Transform/Reductor/Energy.hpp"
#include "QuICC/Transform/Reductor/EnergyD1.hpp"

namespace QuICC {

namespace Transform {

   void DefaultCartesianChebyshevMap::operator()(MapType& m) const
   {
#ifdef QUICC_HAS_CUDA_BACKEND
      //using backend_t = Fft::Chebyshev::LinearMap::viewGpu_t;
      using backend_t = Fft::Chebyshev::LinearMap::base_t;
#else
      //using backend_t = Fft::Chebyshev::LinearMap::viewCpu_t;
      using backend_t = Fft::Chebyshev::LinearMap::base_t;
#endif

      // Create projectors
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::D1<backend_t>>(m, Backward::D1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::D2<backend_t>>(m, Backward::D2::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::P<backend_t>>(m, Backward::P::id());

      // Create integrators
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::P<backend_t>>(m, Forward::P::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I2<backend_t>>(m, Forward::I2P::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I4<backend_t>>(m, Forward::I4P::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I2D1<backend_t>>(m, Forward::I2D1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I4D1<backend_t>>(m, Forward::I4D1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I2_I2D1<backend_t>>(m, Forward::I2ZI2D1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::I4D1_I2<backend_t>>(m, Forward::I4D1ZI2::id());

      // Create reductors
      this->addOperator<Fft::Chebyshev::LinearMap::Reductor::Energy<backend_t>>(m, Reductor::Energy::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Reductor::EnergyD1<backend_t>>(m, Reductor::EnergyD1::id());
   }
 
} // Transform
} // QuICC
