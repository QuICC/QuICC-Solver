/**
 * @file Mpi.hpp
 * @brief Register MPI initialize and finalize
 */

#ifndef QUICC_ENVIRONMENT_REGISTER_MPI_HPP
#define QUICC_ENVIRONMENT_REGISTER_MPI_HPP

#ifndef QUICC_MPI

namespace QuICC {

namespace Environment {

namespace Register {

   /**
    * @brief MPI is not enabled
    */
   struct Mpi
   {
      /**
       * @brief mpi is not enabled
       */
      static const bool isActive = false;

      /**
       * brief ID in initializer map
       */
      static const int initId = -9999;

      /**
       * brief ID in finalizer map
       */
      static const int finalId = -9999;
   };
}
}
}

#else

// System includes
//
#include <mpi.h>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Environment {

namespace Register {

   /**
    * @brief Register MPI initialize + finalize
    */
   class Mpi
   {
      public:
         /**
          * @brief mpi is not enabled
          */
         static const bool isActive = true;

         /**
          * brief Priority in initializer map
          */
         static const int priority;

         /**
          * brief ID in initializer map
          */
         static const int initId;

         /**
          * brief ID in finalizer map
          */
         static const int finalId;

         /**
          * brief Initializer
          */
         static void initializer();

         /**
          * brief Initializer
          */
         static void finalizer();

      protected:

      private:
         /**
          * @brief Constructor
          */
         Mpi();

         /**
          * @brief Destructor
          */
         ~Mpi();
   };

}
}
}

#endif // QUICC_ENVIRONMENT_REGISTER_MPI_HPP
#endif // QUICC_MPI
