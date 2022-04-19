/**
 * @file Hdf5.hpp
 * @brief Register HDF5 initialize and finalize
 */

#ifndef QUICC_ENVIRONMENT_REGISTER_HDF5_HPP
#define QUICC_ENVIRONMENT_REGISTER_HDF5_HPP

#ifndef QUICC_LARGEIO_HDF5

namespace QuICC {

namespace Environment {

namespace Register {

   /**
    * @brief Hdf5 is not enabled
    */
   struct Hdf5
   {
      /**
       * @brief HDF5 is not enabled
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
#include <hdf5.h>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Environment {

namespace Register {

   /**
    * @brief Register HDF5 initialize + finalize
    */
   class Hdf5
   {
      public:
         /**
          * @brief Environment is active
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
         Hdf5();

         /**
          * @brief Destructor
          */
         ~Hdf5();
   };

}
}
}

#endif // QUICC_LARGEIO_HDF5
#endif // QUICC_ENVIRONMENT_REGISTER_HDF5_HPP
