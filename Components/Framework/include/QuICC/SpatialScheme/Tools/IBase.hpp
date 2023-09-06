/** 
 * @file IBase.hpp
 * @brief Implementation of the base tools for the spatial schemes
 */

#ifndef QUICC_SPATIALSCHEME_TOOLS_IBASE_HPP
#define QUICC_SPATIALSCHEME_TOOLS_IBASE_HPP

// System includes
//
#include <vector>
#include <map>

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/Splitting.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   /**
    * @brief Implementation of the base tools for the spatial schemes
    */
   class IBase
   {
      public:
         /**
          * @brief Default ctor
          */
         IBase() = default;

         /**
          * @brief Default dtor
          */
         ~IBase() = default;

         /**
          * @brief Compute backward truncation
          */
         virtual int truncationBwd(const int nN, const int k) = 0;

         /**
          * @brief Compute forward truncation
          */
         virtual int truncationFwd(const int nN, const int k) = 0;

         /**
          * @brief Compute index
          */
         virtual int index(const int nN, const int k) = 0;

         /**
          * @brief Compute total number of modes
          */
         virtual int totalModes(const int n2D, const int n3D, const Splitting::Locations::Id split) = 0;

         /**
          * @brief Check if chosen resolution is optimal
          */
         virtual bool isOptimal(const int nN, const int maxL) = 0;

         /**
          * @brief Build map of balanced indexes for a generic spatial schemes
          */
         virtual void buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>& n0, const std::vector<int>& nN) = 0;

         /**
          * @brief Fill the indexes for 2D and 3D for regular spatial schemes
          */
         void fillIndexes2D3D(std::vector<ArrayI>& idx2D, ArrayI& idx3D, const std::multimap<int,int>& modes);

         /**
          * @brief Fill the indexes for 1D for regular spatial schemes
          */
         void fillIndexes1D(std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, const ArrayI& idx3D, const int nF1D, const int nB1D);
      protected:
         /**
          * @brief Build map of indexes for a generic spatial schemes
          */
         void buildMap(std::multimap<int,int>& modes, const int i0, const int iN, const std::vector<int>& j0, const std::vector<int>& jN, const int c0, const int cN);

         /**
          * @brief Distribute ML load
          */
         void binMLLoad(std::vector<int>& binOrders, const int nL, const int nM, const int id, const int bins);

         /**
          * @brief Compute a simple balanced split of the elements with regard to the given number of parts
          *
          * @param n0         Output start index
          * @param nN         Output number of indexes
          * @param tot        Total number of indexes
          * @param parts      Number of parts to split total into
          * @param id         ID of the CPU
          * @param allowEmpty Allow some parts to be empty
          */
         void balancedSplit(int &n0, int &nN, const int tot, const int parts, const int id, const bool allowEmpty = false);

   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_IBASE_HPP
