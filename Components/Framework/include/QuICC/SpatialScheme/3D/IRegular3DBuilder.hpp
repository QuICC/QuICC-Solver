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
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/Resolutions/Resolution.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"

namespace QuICC {

namespace SpatialScheme {

   /**
    * @brief Implementation of a generic regular 3D scheme
    */
   class IRegular3DBuilder: public IBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim Dimension truncations 
          */
         explicit IRegular3DBuilder(const ArrayI& dim, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~IRegular3DBuilder() = default;

         /**
          * @brief Get spatial scheme dimensions
          */
         virtual ArrayI resolution() const;

         /**
          * @brief Create indexes for a possibly restricted set
          *
          * @param transId Transform ID
          * @param fwd1D   Storage for forward indexes of first dimension
          * @param bwd1D   Storage for backward indexes of first dimension
          * @param idx2D   Storage for the indexes of second dimension
          * @param idx3D   Storage for forward indexes of third dimension
          * @param id      ID of the bin
          * @param bins    Total number of bins (useful to build efficient pairs)
          * @param n0      Starting index of restricted set
          * @param nN      Length of restricted set
          * @param flag    Flag to specify location of splitting
          */
         virtual int fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id = ArrayI(), const ArrayI& bins = ArrayI(), const ArrayI& n0 = ArrayI(), const ArrayI& nN = ArrayI(), const Splitting::Locations::Id flag = Splitting::Locations::NONE);

         /**
          * @brief Get total of splittable indexes 
          *
          * @param transId Transform ID
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag);
         
      protected:
         /**
          * @brief Compute mode distribution for serial algorithm
          *
          * @param modes   Map of all modes
          */
         void splitSerial(std::multimap<int,int>& modes, const Dimensions::Transform::Id transId);

         /**
          * @brief Compute mode distribution for single splitting on 1D algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitSingle1D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

         /**
          * @brief Compute mode distribution for single splitting on 2D algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitSingle2D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

         /**
          * @brief Compute mode distribution for coupled 2D algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitCoupled2D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

         /**
          * @brief Compute mode distribution for tubular 2D decomposition algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitTubular(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

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
