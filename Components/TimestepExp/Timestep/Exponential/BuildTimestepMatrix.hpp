/**
 * @file BuildTimestepMatrixWrapper.hpp
 * @brief Implementation of functors used in Interface
 */

#pragma once

// System includes
//
#include <memory>
#include <map>

// Project includes
//
#include "QuICC/Equations/CouplingInformation.hpp"
#include "QuICC/NonDimensional/INumber.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

/**
 * @brief Wrapper to build timestepping matrices
 */
void buildTimestepMatrix(std::map<std::size_t, DecoupledZSparse>& ops, FieldComponents::Spectral::Id comp, const int idx, const Resolution& res, const Model::IModelBackend& backend, const Equations::CouplingInformation cinfo, const std::map<std::size_t, std::size_t>& bcIds, const std::map<std::size_t, NonDimensional::SharedINumber>& eqParams);

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
