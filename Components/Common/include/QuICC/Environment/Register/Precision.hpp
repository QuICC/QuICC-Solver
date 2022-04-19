/**
 * @file Precision.hpp
 * @brief Register multiple precision initialize and finalize
 */

#ifndef QUICC_ENVIRONMENT_REGISTER_PRECISION_HPP
#define QUICC_ENVIRONMENT_REGISTER_PRECISION_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

namespace Environment {

namespace Register {

   /**
    * @brief Register multiple precision initialize + finalize
    */
   class Precision
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
         Precision();

         /**
          * @brief Destructor
          */
         ~Precision();
   };

}
}
}

#endif // QUICC_ENVIRONMENT_REGISTER_PRECISION_HPP
