/**
 * @file ArnoldiKrylov.hpp
 * @brief Arnoldi iteration for building Krylov subspace
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_ARNOLDIKRYLOV_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_ARNOLDIKRYLOV_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

   /**
    * @brief Arnoldi iteration for building Krylov subspace
    */
   class ArnoldiKrylov
   {
      /**
       * @brief ctor
       */
      ArnoldiKrylov();
      
      /**
       * @brief dtor
       */
      virtual ~ArnoldiKrylov() = default;
   };

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_ARNOLDIKRYLOV_HPP
