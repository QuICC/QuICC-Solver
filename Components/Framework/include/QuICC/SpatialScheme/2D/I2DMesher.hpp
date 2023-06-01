/**
 * @file I2DMesher.hpp
 * @brief Implementation of the 2D spatial scheme mesher
 */

#ifndef QUICC_SPATIALSCHEME_I2DMESHER_HPP
#define QUICC_SPATIALSCHEME_I2DMESHER_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/SpatialScheme/IMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the 2D spatial scheme mesher
    */
   class I2DMesher: public IMesher
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dims Dimension of the domain
          */
         explicit I2DMesher(const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~I2DMesher() = default;

         /**
          * @brief Size of physical grid 3D
          */
         int nPhys3D() const final;

         /**
          * @brief Size of spectral expansion 3D
          */
         int nSpec3D() const final;

         /**
          * @brief Size of dealiased spectral expansion 3D
          */
         int nDealias3D() const final;

      protected:

      private:
   };
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_I2DMESHER_HPP
