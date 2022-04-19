/**
 * @file Hasher.hpp
 * @brief Simple interface to STL hashing function
 */

#ifndef QUICC_HASHER_HPP
#define QUICC_HASHER_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

   /**
    * @brief Simple interface to STL hashing function
    */
   class Hasher
   {
      public:
         /**
          * @brief Compute unique id
          */
         static std::size_t makeId(const std::string tag);
      
      protected:

      private:
         /**
          * @brief Constructor
          */
         Hasher();

         /**
          * @brief Destructor
          */
         virtual ~Hasher(); 
   };

}

#endif // QUICC_HASHER_HPP
