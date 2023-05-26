/** 
 * @file Uniform3D1D.cpp
 * @brief Source of the tools for uniform + uniform + uniform spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/Enums/Splitting.hpp"
#include "QuICC/SpatialScheme/Tools/Uniform3D1D.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int Uniform3D1D::truncationFwd(const int nN, const int l)
   {
      return nN;
   }

   int Uniform3D1D::truncationBwd(const int nN, const int l)
   {
      return nN;
   }

   int Uniform3D1D::index(const int i, const int k)
   {
      return i;
   }

   bool Uniform3D1D::isOptimal(const int nN, const int maxL)
   {
      return true;
   }

   int Uniform3D1D::totalModes(const int n2D, const int n3D, Splitting::Locations::Id split)
   {
      int nModes = 0;
      if(split == Splitting::Locations::FIRST)
      {
         nModes = n2D*n3D;
      }
      else if(split == Splitting::Locations::SECOND)
      {
         nModes = n2D;
      }
      else if(split == Splitting::Locations::BOTH)
      {
         nModes = n3D;
      }
      else
      {
         throw std::logic_error("Splitting location is not implemented");
      }

      return nModes;
   }

   void Uniform3D1D::buildBalancedMap(std::multimap<int,int>& modes, const int n2D, const int n3D, const std::vector<int>& id, const std::vector<int>& bins, const std::vector<int>&, const std::vector<int>&)
   {
      // Create start indexes and length
      std::vector<int> i0(2, 0);
      std::vector<int> iN(2, 0);

      // Get size of the splittable dimension(s)
      int tot = this->totalModes(n2D, n3D, Splitting::Locations::SECOND);

      // Build a simple balanced split
      this->balancedSplit(i0.at(0), iN.at(0), tot, bins.at(1), id.at(1), false);

      // Get size of the splittable dimension(s)
      tot = this->totalModes(n2D, n3D, Splitting::Locations::BOTH);
      tot *= iN.at(0);

      // Compute a balanced splitting
      this->balancedSplit(i0.at(1), iN.at(1), tot, bins.at(0), id.at(0), false);

      if(iN.at(0) > 0 || iN.at(1) > 0)
      {
         int& p0 = i0.at(0);
         int& pN = iN.at(0);
         int& r0 = i0.at(1);
         int& rN = iN.at(1);
         //
         // The compatible dimension is the second dimension, we will need to reorder the indexes
         //

         // Get bottom row offset
         int b0 = r0 % pN;

         // Get the top row length
         int tN = (rN-(pN-b0)) % pN;
         // WARNING:  tN = 0 means that it is a full row!
         if(tN == 0)
         {
            tN = pN;
         }

         // Compute shifted offset
         int s0 = r0/pN;
         // ... and shifted size (+1 because of bottom row)
         int sN = static_cast<int>(std::ceil(static_cast<double>(rN-(pN-b0))/static_cast<double>(pN))) + 1;

         // Create start indexes and length
         std::vector<int> n0(sN+1,p0);
         std::vector<int> nN(sN+1,pN);
         n0.at(0) = s0;
         nN.at(0) = sN;

         // Special treatment for bottom row
         n0.at(1) = (p0 + b0);
         nN.at(1) = pN - b0;

         // Special treatment for top row
         n0.back() = p0;
         nN.back() = tN;

         int k0 = n0.front();
         int kN = nN.front();
         std::vector j0(n0.end()-kN, n0.end());
         std::vector jN(nN.end()-kN, nN.end());
         int c0 = 0;
         int cN = n2D*n3D;

         this->buildMap(modes, k0, kN, j0, jN, c0, cN);
      }
   }

} // Tools
} // SpatialScheme
} // QuICC
