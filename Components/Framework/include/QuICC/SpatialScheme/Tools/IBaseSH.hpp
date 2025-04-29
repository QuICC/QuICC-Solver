/**
 * @file IBaseSH.hpp
 * @brief Implementation of the tools for the generic + spherical harmonics schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_IBASESH_HPP
#define QUICC_SPATIALSCHEME_TOOLS_IBASESH_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/SpatialScheme/Tools/IBase.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the generic + spherical harmonic spatial schemes
    */
   class IBaseSH: public IBase
   {
      public:
         /**
          * @brief Default ctor
          */
         IBaseSH() = default;

         /**
          * @brief ctor with explicit min truncation
          *
          * @param min  Minimal truncation
          */
         IBaseSH(const int min);

         /**
          * @brief Default dtor
          */
         ~IBaseSH() = default;

         /**
          * @brief Compute total number of modes
          */
         int totalModes(const int n2D, const int n3D, const Splitting::Locations::Id split) final;

         /**
          * @brief Build load balance map of indexes for a generic spherical harmonic spatial schemes
          */
         virtual void buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>& n0, const std::vector<int>& nN) override;

      protected:
         /**
          * @brief Compute number of harmonic orders for given harmonic degree l
          *
          * Output value is bounded by max number of harmonic orders
          */
         int nM(const int l, const int nM);

         /**
          * @brief Compute number of harmonic degrees for given harmonic order m
          *
          * Output value is bounded by max number of harmonic degrees
          */
         int nL(const int m, const int nL);

         /**
          * @brief Compute the total number of spherical harmonics for given max degrees and orders
          */
         int nHarmonics(const int nL, const int nM);

         /**
          * @brief Original l map balancing
          */
         void balanceOriginal(std::multimap<int,int>& modes, std::multimap<int,int>& lmap, const int id, std::vector<int>& binLoad);
         void balanceAlternate(std::multimap<int,int>& modes, std::multimap<int,int>& lmap, const int id, std::vector<int>& binLoad);
         void balance(std::multimap<int,int>& modes, std::multimap<int,int>& lmap, const int id, std::vector<int>& binLoad, std::list<int>& keys, const int kind);
         std::list<int> alternatingUpDown(std::multimap<int,int>& lmap, const std::size_t bins);
         std::list<int> oldOrdering(std::multimap<int,int>& lmap, const std::size_t bins);
   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_IBASESH_HPP
