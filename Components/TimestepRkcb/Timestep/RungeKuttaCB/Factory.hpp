/**
 * @file Factory.hpp
 * @brief Factor for Runge-Kutta CB schemes
 */

#ifndef QUICC_TIMESTEP_RUNGEKUTTACB_FACTORY_HPP
#define QUICC_TIMESTEP_RUNGEKUTTACB_FACTORY_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Timestep/Interface.hpp"

namespace QuICC {

namespace Timestep {

/// @brief This namespace provides Runge-Kutta timestepping schemes based on CB15
namespace RungeKuttaCB {

std::shared_ptr<Timestep::Interface> makeInterface(const std::size_t schemeId,
   const MHDFloat time, const Matrix& cfl, const MHDFloat maxError,
   const Timestep::Interface::ScalarEquation_range& scalEq,
   const Timestep::Interface::VectorEquation_range& vectEq,
   Pseudospectral::Coordinator& pseudo);

}
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_RUNGEKUTTACB_FACTORY_HPP
