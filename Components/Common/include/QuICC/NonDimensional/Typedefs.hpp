/**
 * @file Typedefs.hpp
 * @brief Some general typedefs for nondimensional numbers
 */

#ifndef QUICC_NONDIMENSIONAL_TYPEDEFS_HPP
#define QUICC_NONDIMENSIONAL_TYPEDEFS_HPP

// Configuration includes
//

// System includes
//
#include <map>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/NonDimensional/INumber.hpp"

namespace QuICC {

namespace NonDimensional {

   /// Typedef for map of nondimensional numbers
   typedef std::map<std::size_t, SharedINumber>   NdMap;

}
}

#endif // QUICC_NONDIMENSIONAL_TYPEDEFS_HPP
