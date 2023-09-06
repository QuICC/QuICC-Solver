/**
 * @file Kokkos.cpp
 * @brief Kokkos library fixture
 */

// System includes
//
#ifdef QUICC_USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

// Project includes
//
#include "Kokkos.hpp"

namespace QuICC {

namespace ExternalLibrary {

   Kokkos::Kokkos()
   {
      #ifdef QUICC_USE_KOKKOS
      ::Kokkos::initialize();
      #endif
   }

   Kokkos::~Kokkos()
   {
      #ifdef QUICC_USE_KOKKOS
      ::Kokkos::finalize();
      #endif
   }

   Kokkos& Kokkos::getInstance()
   {
      static Kokkos instance;
      return instance;
   }


} // namespace QuICC
} // namespace ExternalLibrary

