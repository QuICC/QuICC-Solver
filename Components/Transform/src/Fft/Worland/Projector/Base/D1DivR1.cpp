/**
 * @file D1DivR1.cpp
 * @brief Source of the implementation of the Worland D1 R^(-1) projector
 * @brief Formula: D1(W_n^l/r) = (norm)^(-1) * r^(l-2) * [ (l-1)*P_n^(alpha,beta) + 2 r^2 (l+n) P_(n-1)^(alpha+1,beta+1) ] 
 */

// External includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Worland/Projector/Base/D1DivR1.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   D1DivR1<base_t>::D1DivR1()
   {
      this->setProfileTag();
   }

   void D1DivR1<base_t>::initBackend() const
   {
      throw std::logic_error(
         "FFT version of D1DivR1 operator is WIP.");
         // specifically, I am not sure that two applications of lowerBeta at the end
         // works properly for l_in = 1

      // l = 0 mode is set to sero: output would be l=-2, non-physical
      std::set<int> filter = {0};
      this->mBackend.setZFilter(filter);

      int lshift = -2;  // The derivative d/dr(1/r) maps Worland-l to Worland-(l-2).
      int extraN = 1;   // The working spectral size is mSpecSize + extraN + lmax/2.
                        // The +1 extra mode is needed because copy(extra, main, -1) shifts
                        // coefficients by one index: extra[n] = main[n+1], so the buffer
                        // must have room for one extra row beyond the normal spectral size.
      bool onlyShiftParity = true;  // If false, lshift would immediately offset all stored l-values by -1.
                                    // The rule is: whenever the intermediate Jacobi operations need the original input l, 
                                    //use onlyShiftParity=true and let lowerBeta apply the decrement(s) explicitly at the end.
      this->mBackend.init(*this->mspSetup, lshift, extraN, onlyShiftParity);  // Reads the setup (list of l-values, spectral size, grid size), partitions
                                                                              // all modes into even/odd DCT buckets using the parity of (l + lshift),
                                                                              // allocates the primary input/output buffers, precomputes the Jacobi
                                                                              // shift matrix J used by all buildShift* calls, and precomputes the
                                                                              // banded matrix pairs for backwardWorland.
      this->mBackend.addStorage(1, 0); // Allocates 1 extra input buffer (the `extra` slot, index 1) and 0
                                       // extra output buffers. This is the scratch space for the second term
                                       // of the derivative. Without this, only slot 0 (main) would exist.
   }

   void D1DivR1<base_t>::computeWorlandExpansion(const bool isEven) const
   {
      const int main = 0;  // index into the backend's buffer array: slot 0 is the primary
                           // working buffer, holds the input spectral coefficients c_n throughout
      const int extra = 1; // slot 1 is the extra buffer allocated by addStorage(1,0)
                           // "1 extra mode" in initBackend means we need exactly one
                           // additional buffer to hold the second term while we build it


      this->mBackend.scaleC(1.0/std::sqrt(Math::PI), isEven);  // multiply all c_n by 1/sqrt(pi). Bridges orthonormal Worland norm
                                                               // to the un-normalized Jacobi representation that backwardWorland expects.

      // SECOND TERM                                                         
      // Copy expansion shifting indexes by -1 and scaling
      this->mBackend.copy(extra, main, -1, isEven);   // Copy main -> extra, shifting the n-index by -1: extra[n] = main[n+1].
                                                      // This implements P_n -> P_{n-1} from the Jacobi derivative formula.
                                                      // After this, extra holds coefficients of P_{n-1}^{(-1/2, l-1/2)}.
      this->mBackend.scaleD(isEven, 0, extra);  // Multiply extra[n] by 2*(l+n+1) * sqrt((n+1)/(n+l+1)).
                                                // This is the derivative coefficient: d/dx P_n^{(a,b)} = (n+a+b+1)/2 * P_{n-1}^{(a+1,b+1)}
                                                // with a=-1/2, b=l-1/2 gives n+l. The sqrt factor converts between
                                                // orthonormal norms of P^{(a,b)} and P^{(a+1,b+1)}.
                                                // After this, extra holds the correctly scaled P_{n-1}^{(1/2, l+1/2)} coefficients.
      this->mBackend.lshift(extra, 1, isEven);  // Pure metadata: tell the backend the current beta of extra is (l+1)-1/2 = l+1/2,
                                                // not l-1/2. Needed so the next two shift operations build their matrices with
                                                // the correct beta. No arithmetic on the coefficients.
      this->mBackend.lowerAlpha(0.5, isEven, extra, 1.0);   // Shift alpha: P^{(+1/2, l+1/2)} -> P^{(-1/2, l+1/2)}.
      this->mBackend.lowerR2Beta(-0.5, isEven, extra, 1.0); // Lower beta by 1 and divide by r^2
                                                            // // r^2 * P^{(-1/2, l+1/2)} -> P^{(-1/2, l-1/2)}

      // FIRST TERM   
      // Scaling original expansion
      this->mBackend.scaleALPY(1.0, -1.0, isEven);  // Multiply coefficients by a*l + y = 1.0*l - 1.0

      // Add second term
      this->mBackend.add(main, extra, 0, isEven);  // main += extra (no n-shift).
                                                   // Both terms are now in the same family P^{(-1/2, l-1/2)} with prefactor r^{l-1}.
                                                   // main = (l-1)*P_n^{(-1/2,l-1/2)} + 2r^2*(n+l)*P_{n-1}^{(1/2,l+1/2)}  [fully reduced]
      this->mBackend.lowerBeta(-0.5, isEven);   // lower beta by 1. First argument is alpha
                                                // Change basis from P^{(-1/2, l-1/2)} to P^{(-1/2, (l-1)-1/2)}.
      this->mBackend.lowerBeta(-0.5, isEven);   // need to do this twice: output l is l_in-2


      this->mBackend.backwardWorland(isEven);   // Evaluates the Worland-(l-1) expansion at physical grid points.
                                                // Output: f'(r_j) at all quadrature nodes.
   }

   void D1DivR1<base_t>::applyPreOperator(const Matrix& in, const bool isEven) const
   {
      this->mBackend.input(in, isEven, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void D1DivR1<base_t>::applyPostOperator(Matrix& rOut, const bool isEven) const
   {
      this->mBackend.output(rOut, isEven);
   }

   void D1DivR1<base_t>::applyPreOperator(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mBackend.input(in, isEven, useReal, true);
      this->computeWorlandExpansion(isEven);
      this->mBackend.io(isEven);
   }

   void D1DivR1<base_t>::applyPostOperator(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mBackend.output(rOut, isEven, useReal);
   }

}
}
}
}
}
