/**
 * file DefaultCylinderWorlandMap.cpp
 * @brief Source of the implementation of the Worland transform in a sphere
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/DefaultCylinderWorlandMap.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/P.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1_P.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1D1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/CylLaplh.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/CylLaplh_DivR1D1R1.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/DivR1CylLaplh_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1CylLaplh.hpp"
#include "QuICC/Transform/Poly/Worland/Projector/D1CylLaplh_D1DivR1D1R1.hpp"

#include "QuICC/Transform/Poly/Worland/Integrator/P.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I4DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I4DivR1D1R1_I2.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I6DivR1_Zero.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I6DivR1D1R1_I4.hpp"
#include "QuICC/Transform/Poly/Worland/Integrator/I6CylLaplh_I4D1R1.hpp"

#include "QuICC/Transform/Forward/P.hpp"
#include "QuICC/Transform/Forward/I4Overr1Pm.hpp"
#include "QuICC/Transform/Forward/I4Overr1D1R1ZI2.hpp"
#include "QuICC/Transform/Forward/I6Overr1Pm.hpp"
#include "QuICC/Transform/Forward/I6Overr1D1R1ZI4.hpp"
#include "QuICC/Transform/Forward/I6LaplhZI4D1R1.hpp"

#include "QuICC/Transform/Backward/P.hpp"
#include "QuICC/Transform/Backward/Overr1Pm.hpp"
#include "QuICC/Transform/Backward/D1.hpp"
#include "QuICC/Transform/Backward/D1ZP.hpp"
#include "QuICC/Transform/Backward/Overr1D1R1.hpp"
#include "QuICC/Transform/Backward/Laplh.hpp"
#include "QuICC/Transform/Backward/LaplhZOverr1D1R1.hpp"
#include "QuICC/Transform/Backward/Overr1LaplhPm.hpp"
#include "QuICC/Transform/Backward/D1Laplh.hpp"
#include "QuICC/Transform/Backward/D1LaplhZD1Overr1D1R1.hpp"

namespace QuICC {

namespace Transform {

   void DefaultCylinderWorlandMap::operator()(MapType& m) const
   {
      using namespace Poly::Worland;
      using backend_t = base_t;

      // Reserve storage for the projectors
      this->addOperator<Projector::P<backend_t>>(m, Backward::P::id());
      this->addOperator<Projector::DivR1_Zero<backend_t>>(m, Backward::Overr1Pm::id());
      this->addOperator<Projector::D1<backend_t>>(m, Backward::D1::id());
      this->addOperator<Projector::D1_P<backend_t>>(m, Backward::D1ZP::id());
      this->addOperator<Projector::DivR1D1R1<backend_t>>(m, Backward::Overr1D1R1::id());
      this->addOperator<Projector::CylLaplh<backend_t>>(m, Backward::Laplh::id());
      this->addOperator<Projector::CylLaplh_DivR1D1R1<backend_t>>(m, Backward::LaplhZOverr1D1R1::id());
      this->addOperator<Projector::DivR1CylLaplh_Zero<backend_t>>(m, Backward::Overr1LaplhPm::id());
      this->addOperator<Projector::D1CylLaplh<backend_t>>(m, Backward::D1Laplh::id());
      this->addOperator<Projector::D1CylLaplh_D1DivR1D1R1<backend_t>>(m, Backward::D1LaplhZD1Overr1D1R1::id());

      // Reserve storage for the weighted projectors
      this->addOperator<Integrator::P<backend_t>>(m, Forward::P::id());
      this->addOperator<Integrator::I4DivR1_Zero<backend_t>>(m, Forward::I4Overr1Pm::id());
      this->addOperator<Integrator::I4DivR1D1R1_I2<backend_t>>(m, Forward::I4Overr1D1R1ZI2::id());
      this->addOperator<Integrator::I6DivR1_Zero<backend_t>>(m, Forward::I6Overr1Pm::id());
      this->addOperator<Integrator::I6DivR1D1R1_I4<backend_t>>(m, Forward::I6Overr1D1R1ZI4::id());
      this->addOperator<Integrator::I6CylLaplh_I4D1R1<backend_t>>(m, Forward::I6LaplhZI4D1R1::id());
   }

}
}
