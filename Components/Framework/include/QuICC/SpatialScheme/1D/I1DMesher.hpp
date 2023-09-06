/**
 * @file I1DMesher.hpp
 * @brief Implementation of the 1D spatial scheme mesher
 */

#ifndef QUICC_SPATIALSCHEME_I1DMESHER_HPP
#define QUICC_SPATIALSCHEME_I1DMESHER_HPP

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
    * @brief Implementation of the 1D spatial scheme mesher
    */
   class I1DMesher: public IMesher
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dims Dimension of the domain
          */
         explicit I1DMesher(const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~I1DMesher() = default;

         /**
          * @brief Size of physical grid 2D
          */
         int nPhys2D() const final;

         /**
          * @brief Size of spectral expansion 2D
          */
         int nSpec2D() const final;

         /**
          * @brief Size of dealiased spectral expansion 2D
          */
         int nDealias2D() const final;

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

#endif // QUICC_SPATIALSCHEME_I1DMESHER_HPP
