/**
 * @file CreateInfoFunctor.hpp
 * @brief 
 */

#pragma once

// System includes
//

// Project includes
//
#include "Timestep/Exponential/TimestepperInfo.hpp"
#include "QuICC/Equations/CouplingInformation.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

namespace Functors {

/**
 * @brief Gather info for Timestepper
 */
TimestepperInfo createInfo(const Equations::CouplingInformation& cinfo, const std::size_t idx, const std::size_t fieldIndex);

inline TimestepperInfo createInfo(const Equations::CouplingInformation& cinfo, const std::size_t idx, const std::size_t fieldIndex)
{
   TimestepperInfo info;
   info.isComplex = false;
   info.fieldIndex = fieldIndex;
   info.solverIndex = 0;
   info.rows = 2*cinfo.systemN(idx)*cinfo.rhsCols(idx);
   info.cols = 1;
   info.blockN = 2*cinfo.galerkinN(idx)*cinfo.rhsCols(idx);


   return info;
}


} // namespace Functors
} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
