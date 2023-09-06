/**
 * @file Constants.hpp
 * @brief Some constants used in timestepping operations
 */

#ifndef QUICC_TIMESTEP_CONSTANTS_HPP
#define QUICC_TIMESTEP_CONSTANTS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Timestep {

   /**
    * @brief Special CFL location value for fixed step condition
    */
   static const MHDFloat FIXEDSTEP_LOCATION = -100;

   /**
    * @brief Special CFL location value for max step condition
    */
   static const MHDFloat MAXSTEP_LOCATION = -101;

   /**
    * @brief Special CFL location value for min step condition
    */
   static const MHDFloat MINSTEP_LOCATION = -102;

   /**
    * @brief Special CFL location value for error step condition
    */
   static const MHDFloat ERROR_LOCATION = -200;

   /**
    * @brief Limit for max timestep size
    */
   static const MHDFloat LIMIT_MAXSTEP = 1e-1;

   /**
    * @brief Limit for min timestep size
    */
   static const MHDFloat LIMIT_MINSTEP = 1e-10;

   /**
    * @brief Only just if timestep increase is more than 5%
    */
   static const MHDFloat FIVE_PC_WINDOW = 1.05;

   /**
    * @brief Maximum stepsize increase in single step
    *
    * Soederlind, DOI: 10.1023/A:1021160023092
    * Usman and Hall, DOI: 10.1016/S0377-0427(97)00239-2
    */
   static const MHDFloat MAX_STEPSIZE_JUMP = 1.602;

}
}

#endif // QUICC_TIMESTEP_CONSTANTS_HPP
