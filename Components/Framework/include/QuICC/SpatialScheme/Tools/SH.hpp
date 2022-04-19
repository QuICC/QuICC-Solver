/** 
 * @file SH.hpp
 * @brief Implementation of the tools for the spherical harmonics based schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_SH_HPP
#define QUICC_SPATIALSCHEME_TOOLS_SH_HPP

// Configuration includes
//

// System includes
//
#include <vector>
#include <deque>
#include <queue>
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the tools for the spherical harmonics based schemes
    */
   class SH
   {
      public:
         /**
          * @brief Compute number of harmonic orders for given harmonic degree l
          *
          * Output value is bounded by max number of harmonic orders
          */
         static int nM(const int l, const int nM);

         /**
          * @brief Compute number of harmonic degrees for given harmonic order m
          *
          * Output value is bounded by max number of harmonic degrees
          */
         static int nL(const int m, const int nL);

         /**
          * @brief Compute the total number of spherical harmonics for given max degrees and orders
          */
         static int nHarmonics(const int nL, const int nM);

         /**
          * @brief Compute map of all harmonics with harmonic order m as key
          */
         static void buildMMap(std::multimap<int,int>& harmonics, const int nL, const int nM);

         /**
          * @brief Compute map of all harmonics with harmonic degree l as key
          */
         static void buildLMap(std::multimap<int,int>& harmonics, const int nL, const int nM);

         /**
          * @brief Compute map of given harmonic indexes sorted to improved load balance with harmonic degree l as key
          */
         static void buildLHSortedMap(std::multimap<int,int>& harmonics, const int nL, const int nM, const int h0, const int nH);

         /**
          * @brief Compute map of given harmonic indexes sorted to improved load balance with harmonic degree m as key
          */
         static void buildMHSortedMap(std::multimap<int,int>& harmonics, const int nL, const int nM, const int h0, const int nH);

         static void binMLLoad(std::vector<int>& binOrders, const int nL, const int nM, const int id, const int bins);
         /**
          * @brief Compute list of loads for harmonic order m distribution, initialise load per bin and compute optimal load
          */
         static void initMLLoad(std::deque<int>& loadList, std::vector<int>& binLoad, std::queue<int>& optimal, const int nL, const int nM, const int bins);

         static void combineMPairs(std::deque<int>& list, std::multimap<int, int>& regular, const int nM, const int bins);

         static void fillMPairsBins(std::deque<int>& list, std::multimap<int, int>& regular, std::vector<int>& load, const int bins);

         static void fillMRestBins(std::deque<int>& list, std::multimap<int, int>& regular, std::vector<int>& load, std::queue<int>& optimal, const int bins);

         static void convertLoadToOrders(std::multimap<int, int>& regular, const int nL);

         /**
          * @brief Fill the indexes for 1D for SH based spatial schemes
          */
         static void fillIndexes1D(std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, const ArrayI& idx3D, const int nF1D, const int nB1D);
   };
}
}
}

#endif // QUICC_SPATIALSCHEME_TOOLS_SH_HPP
