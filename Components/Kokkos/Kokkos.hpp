/**
 * @file Kokkos.hpp
 * @brief Kokkos library fixture
 */

#ifndef QUICC_KOKKOS_HPP
#define QUICC_KOKKOS_HPP

// System includes
//

// Project includes
//

namespace QuICC {

namespace ExternalLibrary {

   /**
    * @brief Kokkos library fixture
    */
   class Kokkos
   {
      public:
         /**
          * @brief Kokkos library fixture ctor
          */
         Kokkos();

         /**
          * @brief Kokkos library fixture dtor
          */
         ~Kokkos();

         /**
          * @brief return singleton
          */
         static Kokkos& getInstance();
   };

} // namespace QuICC
} // namespace ExternalLibrary

#endif // QUICC_KOKKOS_HPP
