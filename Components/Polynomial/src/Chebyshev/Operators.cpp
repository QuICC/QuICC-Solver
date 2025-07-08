/**
 * @file Tools.cpp
 * @brief Source of the tools for Chebyshev polynomial implementation
 */

// System includes
//
#include <cmath>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Chebyshev/Operators.hpp"

// Project includes
//
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace Polynomial {

namespace Chebyshev {

namespace Operators {

   // Functions needed for the computation of Bmat

   // normalization factor for the Chebyshev polynomials
   static int fT(int n)
   {
      if (n==0) { return 1; }
      else { return 2; }
   }

   // factor to half the first element of the sum in (2) of Saibaba, 2021
   static MHDFloat fc(int n)
   {
      if (n==0) { return 0.5; }
      else { return 1; }
   }

   //Integral of Chebyshev polynomial, Tn over -1<=x<=1
   static MHDFloat IntT(int n)
   {
      if (n==1) { return 0; }
      else { return (std::pow(-1, n) +1)/(1-n*n); }
   }

   // Integral of Tn Tj over -1<=x<=1
   MHDFloat IntTnTj(int n, int j)
   {
      return 0.5*(IntT(n+j)+IntT(Internal::Math::abs(n-j)));
   }

   // ratio  IntTnTj(int n, int j-2) / IntTnTj(int n, int j)
   static MHDFloat ratioIntTnTj2(int n, int j)
   {
      if (IntTnTj(n, j)==0.0) 
      {
         throw std::logic_error("Integral IntTnTj(n, j) is zero!");
      }
      return IntTnTj(n, j-2)/IntTnTj(n, j);
   }

   // ratio of cBar(n,j-2)/cBar(n,j)
   static MHDFloat ratioC2(int n, int j)
   {
      return (fc(j-2)/fc(j) ) * (n+j)/(n-j+2);
   }

   // ratio for the recursive relation to calculate tmpSum in Bmat
   static MHDFloat ratioS(int p, int k, int n0, int n, int alpha)
   {
      return ratioC2(p-k,n0+2*alpha) * ratioIntTnTj2(n,n0+2*alpha);
   }

   // Function to compute the integral of r^p T_n(x)
   void integrateRpTn(Internal::Matrix& iop, const int p, const int nN, const MHDFloat ro, const MHDFloat ri)
   {

      // Chebyshev grid: x = (r-a)/b
      auto a = 0.5*(ro-ri);
      auto b = 0.5*(ro+ri);

      Internal::Array Avec(p,1);
      Internal::Matrix Bmat(p, nN);

      // Avec entries
      Avec(0,0) = std::pow(a, 1 + p);
      for(int k = 0; k<p; k++)
      {
         Avec(k+1,0) = Avec(k) * (b/a) * (p-k)/(k+1);
      }
      
      // Bmat entries
      Bmat.setZero();

      for(int k = 0; k<=p; k++)
      {
         int n0=1;
         if ( (p-k)%2==0 )
         {
            n0=0;
         }
         for(int in = 0; in<(nN-n0); in++)
         {
            int n = n0+2*in;
            auto alphaMax = 0.5*(p - k - n0);
            MHDFloat SalphaMax = std::pow(2, 1-p+k) * fc(n0+2*alphaMax) * IntTnTj(n,n0+2*alphaMax);

            MHDFloat oldTerm = SalphaMax;
            MHDFloat tmpSum = SalphaMax;
            MHDFloat newTerm;

            for(auto alpha = alphaMax; alpha>=1.0; alpha--)
            {
               newTerm = oldTerm * ratioS(p, k, n0, n, alpha);
               tmpSum = tmpSum + newTerm;
               oldTerm = newTerm;
            }
            // Final Matrix
            Bmat(k,n) = fT(n) * tmpSum;
         }
      }

      iop = (Avec.transpose()*Bmat).cast<MHDFloat>().transpose();
   }

}
}
}
}
