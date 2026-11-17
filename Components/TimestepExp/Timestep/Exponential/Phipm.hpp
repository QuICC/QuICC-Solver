/**
 * @file Phipm.hpp
 * @brief phipm algorithm for exponential integrators
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_PHIPM_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_PHIPM_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

   /**
    * @brief phipm algorithm for exponential integrators
    */
   class Phipm
   {
      /**
       * @brief ctor
       */
      Phipm();
      
      /**
       * @brief dtor
       */
      virtual ~Phipm() = default;
   };

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_PHIPM_HPP
