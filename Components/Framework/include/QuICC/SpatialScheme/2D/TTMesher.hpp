/** 
 * @file TTMesher.hpp
 * @brief Implementation of the TT spatial scheme mesher
 */

#ifndef QUICC_SPATIALSCHEME_2D_TTMESHER_HPP
#define QUICC_SPATIALSCHEME_2D_TTMESHER_HPP

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/2D/I2DMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the TT spatial scheme mesher
    */
   class TTMesher: public I2DMesher
   {
      public:
         /**
          * @brief Constructor
          *
          * @param purpose Purpose of the grid
          */
         explicit TTMesher(const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         ~TTMesher() = default;

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
         
      protected:
         /**
          * @brief X grid size
          */
         int mNx;

         /**
          * @brief Y grid size
          */
         int mNy;

      private:
   };
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_2D_TTMESHER_HPP
