/**
 * @file Kiops.hpp
 * @brief KIOPS algorithm for exponential integrators
 */

#ifndef QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP
#define QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace Timestep {

/// @brief This namespace contains exponential timestepping schemes
namespace Exponential {

   /**
    * @brief KIOPS algorithm for exponential integrators
    */
   class Kiops
   {
      public:
         /**
          * @brief ctor
          */
         Kiops();
         
         /**
          * @brief dtor
          */
         virtual ~Kiops() = default;
      private:
         /**
          * @brief Tolerance
          */
         const double mcTol;

         /**
          * @brief Min size Krylov subspace
          */
         const double mcMMin;

         /**
          * @brief Max size of Krylov subspace
          */
         const double mcMMax;

         /**
          * @brief Size of Krylov subspace
          */
         double mM;
   };

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC

#endif // QUICC_TIMESTEP_EXPONENTIAL_KIOPS_HPP
