/**
 * @file Utils.cpp
 * @brief Pseudospectral utils
 */

// System includes
//
#include <stdexcept>
#include <type_traits>
#include <boost/functional/hash.hpp>

// Project includes
//
#include "QuICC/Pseudospectral/Utils.hpp"

namespace QuICC {

namespace Pseudospectral {

namespace details
{
    std::size_t hash_combine(const std::size_t a, const std::size_t b)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, a);
        boost::hash_combine(seed, b);
        return seed;
    }

} // namespace details
} // Pseudospectral
} // QuICC
