/**
 * file DefaultAnnulusChebyshevMap.cpp
 * @brief Source of the implementation of the Chebyshev transform in a annulus
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/DefaultAnnulusChebyshevMap.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY2.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
//NOt IMPLEMENTED YET #include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1D1.hpp"
//NOt IMPLEMENTED YET#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1DivY1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
                                              
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"

#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Energy.hpp"

#include "QuICC/Transform/Forward/P.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1.hpp"
#include "QuICC/Transform/Backward/Overr2.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
//NOt IMPLEMENTED YET#include "QuICC/Transform/Backward/Overr1D1.hpp"
//NOt IMPLEMENTED YET#include "QuICC/Transform/Backward/D1Overr1.hpp"

#include "QuICC/Transform/Reductor/Energy.hpp"

namespace QuICC {

namespace Transform {

   void DefaultAnnulusChebyshevMap::operator()(MapType& m) const
   {
      // Create projectors
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::P>(m, Backward::P::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY1>(m, Backward::Overr1::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::DivY2>(m, Backward::Overr2::id());
      this->addOperator<Fft::Chebyshev::LinearMap::Projector::D1>(m, Backward::D1::id());
      //NOT IMPLEMENTED YET this->addOperator<Fft::Chebyshev::LinearMap::Projector::D1>(m, Backward::Overr1D1::id());
      //NOT IMPLEMENTED YET this->addOperator<Fft::Chebyshev::LinearMap::Projector::D1>(m, Backward::D1Overr1::id());

      // Create integrators
      this->addOperator<Fft::Chebyshev::LinearMap::Integrator::P>(m, Forward::P::id());

      // Create reductors
      this->addOperator<Fft::Chebyshev::LinearMap::Reductor::Energy>(m, Reductor::Energy::id());
   }
 
} // Transform
} // QuICC
