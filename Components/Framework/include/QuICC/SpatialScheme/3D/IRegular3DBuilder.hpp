/**
 * @file IRegular3DBuilder.hpp
 * @brief Implementation of a generic regular 3D scheme
 */

#ifndef QUICC_SPATIALSCHEME_IREGULAR3DBUILDER_HPP
#define QUICC_SPATIALSCHEME_IREGULAR3DBUILDER_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/3D/I3DBuilder.hpp"
#include "QuICC/SpatialScheme/Tools/IBase.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of a generic regular 3D scheme
    */
   class IRegular3DBuilder: public I3DBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim Dimension truncations
          * @param purpose Setup purpose: simulation, visualization
          * @param options Options for builder
          */
         explicit IRegular3DBuilder(const ArrayI& dim, const GridPurpose::Id purpose, const std::map<std::size_t,std::vector<std::size_t>>& options);

         /**
          * @brief Destructor
          */
         virtual ~IRegular3DBuilder() = default;

         /**
          * @brief Get spatial scheme dimensions
          */
         virtual ArrayI resolution() const override;

      protected:
         /**
          * @brief Get truncation tool
          */
         virtual std::shared_ptr<Tools::IBase> truncationTools(const Dimensions::Transform::Id transId) const override;

         /**
          * @brief Compute mode distribution for coupled 2D algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitCoupled2D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

         /**
          * @brief First regular truncation
          */
         int   mI;

         /**
          * @brief Second regular truncation
          */
         int   mJ;

         /**
          * @brief third regular truncation
          */
         int   mK;

      private:
   };

} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_IREGULAR3DBUILDER_HPP
