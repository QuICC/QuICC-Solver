/**
 * @file CreatorImpl.hpp
 * @brief Spatial scheme creator implementation
 */

#ifndef QUICC_SPATIALSCHEME_CREATORIMPL_HPP
#define QUICC_SPATIALSCHEME_CREATORIMPL_HPP

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SpatialScheme/ICreator.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Spatial scheme creator implementation
    */
   template <typename TName> class CreatorImpl: public ICreator
   {
      public:
         /**
          * @brief Constructor
          */
         CreatorImpl();

         /**
          * @brief Destructor
          */
         virtual ~CreatorImpl();

         /**
          * @brief Create spatial scheme shared pointer
          */
         virtual SharedISpatialScheme create(const VectorFormulation::Id formulation, const GridPurpose::Id purpose) const;
   };

   template <typename TName> CreatorImpl<TName>::CreatorImpl()
   {
   }

   template <typename TName> CreatorImpl<TName>::~CreatorImpl()
   {
   }

   template <typename TName> SharedISpatialScheme CreatorImpl<TName>::create(const VectorFormulation::Id formulation, const GridPurpose::Id purpose) const
   {
      std::shared_ptr<TName> spName = std::make_shared<TName>(formulation, purpose);

      return spName;
   }

}
}

#endif // QUICC_SPATIALSCHEME_CREATORIMPL_HPP
