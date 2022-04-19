/**
 * @file Coordinator.hpp
 * @brief Spatial scheme coordinator
 */

#ifndef QUICC_SPATIALSCHEME_COORDINATOR_HPP
#define QUICC_SPATIALSCHEME_COORDINATOR_HPP

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Hasher.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SpatialScheme/ICreator.hpp"
#include "QuICC/SpatialScheme/CreatorImpl.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Spatial scheme coordinator
    */
   class Coordinator
   {
      public:
         /**
          * @brief Map indexes to creators
          */
         static std::map<std::size_t, SharedICreator>&  map();
      
      protected:

      private:
         /**
          * @brief Constructor
          */
         Coordinator();

         /**
          * @brief Destructor
          */
         virtual ~Coordinator(); 
   };
      
   /**
    * @brief Register new spatial scheme
    */
   template <typename TScheme> static std::size_t registerId(const std::string tag);

   template <typename TScheme> std::size_t registerId(const std::string s)
   {
      std::size_t id = Hasher::makeId(s);
      std::shared_ptr<CreatorImpl<TScheme> >  spFactory = std::make_shared<CreatorImpl<TScheme> >();

      Coordinator::map().insert(std::make_pair(id,spFactory));

      return id;
   }

}
}

#endif // QUICC_SPATIALSCHEME_COORDINATOR_HPP
