/** 
 * @file IRegularSHmlBuilder.hpp
 * @brief Implementation of the Regular basis + Spherical Harmonics scheme with m spectral ordering and l spectral transform ordering
 */

#ifndef QUICC_SPATIALSCHEME_IREGULARSHMLBUILDER_HPP
#define QUICC_SPATIALSCHEME_IREGULARSHMLBUILDER_HPP

// Configuration includes
//

// System includes
//

// External includes
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
    * @brief Implementation of Regular basis + Spherical Harmonics scheme with m spectral ordering and l spectral transform ordering
    *
    */
   class IRegularSHmlBuilder: public IBuilder
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim 0: radial, 1: latitudinal, 2: longitudinal
          */
         explicit IRegularSHmlBuilder(const ArrayI& dim, const GridPurpose::Id purpose);

         /**
          * @brief Destructor
          */
         virtual ~IRegularSHmlBuilder();

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

         /**
          * @brief Add index counter to shared resolution
          */
         virtual void addIndexCounter(SharedResolution spRes);
         
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
         void splitSingle2D(std::multimap<int,int>& modes, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

         /**
          * @brief Compute mode distribution for tubular 2D decomposition algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitTubular(std::multimap<int,int>& modes, const ArrayI& id, const ArrayI& bins, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

         /**
          * @brief Regular truncation
          */
         int   mI;

         /**
          * @brief Spherical harmonic degree
          */
         int   mL;

         /**
          * @brief Spherical harmonic order
          */
         int   mM;

      private:
         /**
          * @brief Build the ML map of the spherical harmonics
          */
         void buildMLMap(std::multimap<int,int>& harmonics, const int id, const int bins);
   };
}
}

#endif // QUICC_SPATIALSCHEME_IREGULARSHMLBUILDER_HPP
