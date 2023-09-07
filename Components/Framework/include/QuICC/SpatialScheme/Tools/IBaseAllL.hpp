/** 
 * @file IBaseAllL.hpp
 * @brief Implementation of the base tools for the radial + spherical harmonics schemes with all harmonic degrees gathered
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_IBASEALLL_HPP
#define QUICC_SPATIALSCHEME_TOOLS_IBASEALLL_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/IBaseSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the radial + spherical harmonic spatial schemes with all harmonic degrees gathered
    */
   class IBaseAllL: public IBaseSH
   {
      public:
         /**
          * @brief Default ctor
          */
         IBaseAllL() = default;

         /**
          * @brief Default dtor
          */
         ~IBaseAllL() = default;

         /**
          * @brief Compute forward truncation (calls truncationBwd)
          *
          * @param nN   Reference truncation
          * @param j    second dimension
          * @param k    third dimension
          */
         int truncationFwd(const int nN, const int j, const int k) final;

         /**
          * @brief Build load balance map of indexes for a generic spherical harmonic spatial schemes
          */
         virtual void buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>& n0, const std::vector<int>& nN) override;

      protected:
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_IBASEALLL_HPP
