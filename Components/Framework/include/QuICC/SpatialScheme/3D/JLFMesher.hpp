/** 
 * @file JLFMesher.hpp
 * @brief Implementation of the JLF spatial scheme mesher
 */

#ifndef QUICC_SPATIALSCHEME_3D_JLFMESHER_HPP
#define QUICC_SPATIALSCHEME_3D_JLFMESHER_HPP

// System includes
//

// Project includes
//
#include "QuICC/SpatialScheme/IMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of the JLF spatial scheme mesher
    */
   class JLFMesher: public IMesher
   {
      public:
         /**
          * @brief Constructor
          *
          * @param purpose Purpose of the grid
          */
         explicit JLFMesher(const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         ~JLFMesher() = default;

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
         /**
          * @brief Dealiased spectral expansion 1D
          */
         int mNdealias;

         /**
          * @brief Radial grid size
          */
         int mNr;

         /**
          * @brief Theta grid size
          */
         int mNt;

         /**
          * @brief Phi grid size
          */
         int mNp;

      private:
   };
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_3D_JLFMESHER_HPP
