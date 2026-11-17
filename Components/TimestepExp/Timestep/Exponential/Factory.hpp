/**
 * @file Factory.hpp
 * @brief Exponential timestepper factory
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_FACTORY_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_FACTORY_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Timestep/Interface.hpp"

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

std::shared_ptr<Timestep::Interface> makeInterface(const std::size_t schemeId,
   const MHDFloat time, const Matrix& cfl, const MHDFloat maxError,
   const Timestep::Interface::ScalarEquation_range& scalEq,
   const Timestep::Interface::VectorEquation_range& vectEq,
   Pseudospectral::Coordinator& pseudo);

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_FACTORY_HPP
