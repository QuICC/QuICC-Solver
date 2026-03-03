/**
 * @file IomKrylov.hpp
 * @brief Incomplete orthogonalization method for building Krylov subspace
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_IOMKRYLOV_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_IOMKRYLOV_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

   /**
    * @brief Incomplete orthogonalization method for building Krylov subspace
    */
   class IomKrylov
   {
      /**
       * @brief ctor
       */
      IomKrylov();
      
      /**
       * @brief dtor
       */
      virtual ~IomKrylov() = default;
   };

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_IOMKRYLOV_HPP
