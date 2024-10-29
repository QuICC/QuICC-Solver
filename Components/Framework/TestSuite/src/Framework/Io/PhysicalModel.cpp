/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a sphere
 * (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/Explicit/PhysicalModel.hpp"
#include "Model/Boussinesq/Sphere/RTC/Explicit/ModelBackend.hpp"
#include "QuICC/Model/PyModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Explicit {

std::string PhysicalModel::PYMODULE()
{
   return "boussinesq.sphere.rtc.explicit.physical_model";
}

void PhysicalModel::init()
{
#ifdef QUICC_MODEL_BOUSSINESQSPHERERTC_EXPLICIT_BACKEND_CPP
   IPhysicalModel<Simulation, StateGenerator, VisualizationGenerator>::init();

   this->mpBackend = std::make_shared<ModelBackend>();
#else
   IPhysicalPyModel<Simulation, StateGenerator, VisualizationGenerator>::init();

   this->mpBackend =
      std::make_shared<PyModelBackend>(this->PYMODULE(), this->PYCLASS());
#endif
}

} // namespace Explicit
} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
