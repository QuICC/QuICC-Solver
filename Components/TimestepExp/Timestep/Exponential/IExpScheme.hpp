/**
 * @file IExpScheme.hpp
 * @brief Interface for exponential timestep scheme
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_IEXPSCHEME_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_IEXPSCHEME_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/Timestep/IScheme.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

/**
 * @brief Interface for expoential scheme 
 */
class IExpScheme: public IScheme
{
public:
   /**
    * @brief Constructor
    */
   IExpScheme() = default;

   /**
    * @brief Destructor
    */
   virtual ~IExpScheme() = default;

protected:

private:
};

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_IEXPSCHEME_HPP
