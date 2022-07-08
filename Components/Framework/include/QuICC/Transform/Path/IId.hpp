/**
 * @file IId.hpp
 * @brief Interface for a generic path ID
 */

#ifndef QUICC_TRANSFORM_PATH_IID_HPP
#define QUICC_TRANSFORM_PATH_IID_HPP

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   /**
    * @brief Interface for a generic path ID
    */
   class IId
   {
      public:
         /**
          * @brief Constructor
          */
         IId(const std::string tag);

         /**
          * @brief Destructor
          */
         virtual ~IId();

         /**
          * @brief Tag
          */
         std::string tag() const;
         
      protected:

      private:
         /**
          * @brief Tag
          */
         const std::string mTag;
   };

   /// Typedef for shared_pointer IId
   typedef std::shared_ptr<IId> SharedIId;

}
}
}

#endif // QUICC_TRANSFORM_PATH_IID_HPP
