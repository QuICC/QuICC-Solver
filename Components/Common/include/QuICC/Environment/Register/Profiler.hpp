/**
 * @file Profiler.hpp
 * @brief Register profiler initialize and finalize
 */

#ifndef QUICC_ENVIRONMENT_REGISTER_PROFILER_HPP
#define QUICC_ENVIRONMENT_REGISTER_PROFILER_HPP

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
    * @brief Register profiler initialize + finalize
    */
   class Profiler
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
         Profiler();

         /**
          * @brief Destructor
          */
         ~Profiler();
   };

}
}
}

#endif // QUICC_ENVIRONMENT_REGISTER_PROFILER_HPP
