/**
 * @file IMesher.hpp
 * @brief Implementation of the basic components of the spatial scheme mesher
 */

#ifndef QUICC_SPATIALSCHEME_IMESHER_HPP
#define QUICC_SPATIALSCHEME_IMESHER_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "QuICC/Enums/GridPurpose.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the basic components of the spatial scheme mesher
    */
   class IMesher
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dims Dimension of the domain
          */
         explicit IMesher(const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~IMesher() = default;

         /**
          * @brief Initialize mesher
          */
         virtual void init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Size of physical grid 1D
          */
         virtual int nPhys1D() const = 0;

         /**
          * @brief Size of spectral expansion 1D
          */
         virtual int nSpec1D() const = 0;

         /**
          * @brief Size of dealiased spectral expansion 1D
          */
         virtual int nDealias1D() const = 0;

         /**
          * @brief Size of physical grid 2D
          */
         virtual int nPhys2D() const = 0;

         /**
          * @brief Size of spectral expansion 2D
          */
         virtual int nSpec2D() const = 0;

         /**
          * @brief Size of dealiased spectral expansion 2D
          */
         virtual int nDealias2D() const = 0;

         /**
          * @brief Size of physical grid 3D
          */
         virtual int nPhys3D() const = 0;

         /**
          * @brief Size of spectral expansion 3D
          */
         virtual int nSpec3D() const = 0;

         /**
          * @brief Size of dealiased spectral expansion 3D
          */
         virtual int nDealias3D() const = 0;

      protected:
         /**
          * @brief Main purpose of the grid values
          */
         const GridPurpose::Id mPurpose;

         /**
          * @brief Dimensions
          */
         std::vector<int>  mDims;

      private:
   };
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_IMESHER_HPP
