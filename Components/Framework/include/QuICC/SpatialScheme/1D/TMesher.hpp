/** 
 * @file TMesher.hpp
 * @brief Implementation of the T spatial scheme mesher
 */

#ifndef QUICC_SPATIALSCHEME_1D_TMESHER_HPP
#define QUICC_SPATIALSCHEME_1D_TMESHER_HPP

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/1D/I1DMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the T spatial scheme mesher
    */
   class TMesher: public I1DMesher
   {
      public:
         /**
          * @brief Constructor
          *
          * @param purpose Purpose of the grid
          */
         explicit TMesher(const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         ~TMesher() = default;

         /**
          * @brief Initialize sizes
          *
          * @param dims    Global dimensions
          * @param options Scheme options
          */
         void init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options) final;

         /**
          * @brief Size of physical grid 1D
          */
         int nPhys1D() const final;

         /**
          * @brief Size of spectral expansion 1D
          */
         int nSpec1D() const final;

         /**
          * @brief Size of dealiased spectral expansion 1D
          */
         int nDealias1D() const final;
         
      protected:
         /**
          * @brief X grid size
          */
         int mNx;

      private:
   };
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_1D_TMESHER_HPP
