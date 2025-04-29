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
#include "Types/Typedefs.hpp"
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
          * @brief Named minimal truncation
          */
         enum class MinimalTruncation
         {
            /// Single mode as minimum
            Single = 1,
            /// A minimum of 4 modes is required
            /// to satisfy boundary condition and stable
            /// timestepping
            Triangular = 4,
         };

         /**
          * @brief Default ctor
          */
         IBase();

         /**
          * @brief ctor with explicit min truncation
          *
          * @param min  Minimal truncation
          */
         IBase(const int min);

         /**
          * @brief Default dtor
          */
         ~IBase() = default;

         /**
          * @brief Compute backward truncation
          *
          * @param nN   reference truncation
          * @param j    index of second dimension
          * @param k    index of third dimension
          */
         virtual int truncationBwd(const int nN, const int j, const int k) = 0;

         /**
          * @brief Compute forward truncation
          *
          * @param nN   reference truncation
          * @param j    index of second dimension
          * @param k    index of third dimension
          */
         virtual int truncationFwd(const int nN, const int j, const int k) = 0;

         /**
          * @Brief Minimal truncation
          */
         virtual int min() const;

         /**
          * @brief Compute index
          *
          * @param nN   reference truncation
          * @param j    index of second dimension
          * @param k    index of third dimension
          */
         virtual int index(const int nN, const int j, const int k) = 0;

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
         void fillIndexes2D3D(std::vector<std::vector<int> >& idx2D, std::vector<int>& idx3D, const std::multimap<int,int>& modes);

         /**
          * @brief Fill the indexes for 1D for regular spatial schemes
          */
         void fillIndexes1D(std::vector<std::vector<std::vector<int> > >& fwd1D, std::vector<std::vector<std::vector<int> > >& bwd1D, const std::vector<std::vector<int> >& idx2D, const std::vector<int>& idx3D, const int nF1D, const int nB1D);
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

      protected:
         /**
          * @brief Minimal truncation
          */
         const int mcMin;

   };

} // Tools
} // SpatialScheme
} // QuICC

#endif // QUICC_SPATIALSCHEME_TOOLS_IBASE_HPP
