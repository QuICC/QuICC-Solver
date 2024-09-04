/**
 * @file Tags.hpp
 * @brief Tag types
 */
#pragma once

// External includes
//

// Project includes
//

namespace QuICC {
namespace Graph {

//
// Tags
//

/// @brief Transform data states.
/// From fully physical (PPP) space to fully spectral/modal (MMM) space
/// \todo change stage order to match old nomenclature (i.e. MMM is stage 0)?
enum class Stage
{
    /// physical space, stage 0
    PPP,
    /// modal space, stage 0
    MPP,
    /// physical space, stage 1
    PPM,
    /// modal space, stage 1
    MPM,
    /// physical space, stage 2
    PMM,
    /// modal space, stage 2
    MMM
};

} // namespace Graph
} // namespace QuICC
