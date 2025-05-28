/** 
 * @file FdWLFMesher.hpp
 * @brief Implementation of the FdWLF spatial scheme mesher
 */

#ifndef QUICC_SPATIALSCHEME_3D_FDWLFMESHER_HPP
#define QUICC_SPATIALSCHEME_3D_FDWLFMESHER_HPP

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/3D/xLFMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the FdWLF spatial scheme mesher
    */
   class FdWLFMesher: public xLFMesher
   {
      public:
         /**
          * @brief Constructor
          *
          * @param purpose Purpose of the grid
          */
         explicit FdWLFMesher(const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         ~FdWLFMesher() = default;

         /**
          * @brief Initialize sizes
          *
          * @param dims    Global dimensions
          * @param options Scheme options
          */
         void init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options) final;
         
      protected:

      private:
   };
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_3D_FDWLFMESHER_HPP
