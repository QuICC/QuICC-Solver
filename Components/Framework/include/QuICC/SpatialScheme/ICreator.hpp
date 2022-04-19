/**
 * @file ICreator.hpp
 * @brief Spatial scheme creator interface
 */

#ifndef QUICC_SPATIALSCHEME_ICREATOR_HPP
#define QUICC_SPATIALSCHEME_ICREATOR_HPP

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Spatial scheme creator interface
    */
   class ICreator
   {
      public:
         /**
          * @brief Constructor
          */
         ICreator();

         /**
          * @brief Destructor
          */
         virtual ~ICreator();

         /**
          * @brief Create spatial scheme shared pointer
          */
         virtual SharedISpatialScheme create(const VectorFormulation::Id formulation, const GridPurpose::Id purpose) const = 0;
   };

   /// Typedef for a shared pointer ICreator
   typedef std::shared_ptr<ICreator> SharedICreator;

}
}

#endif // QUICC_SPATIALSCHEME_ICREATOR_HPP
