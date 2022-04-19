/**
 * @file GpuDctTools.cu
 * @brief Source of the Chebyshev GPU utility kernels
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/GpuDctTools.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

namespace Gpu {

   __global__ void buildDCTInput(cufftDoubleReal* out, const cufftDoubleReal *in, const int size, const int blockSize)
   {
      for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
      {
         int c = j*size;
         int k = 2*threadIdx.y;
         for(int i = threadIdx.y; i < size/2; i+=blockDim.y)
         {
            out[i+c] = in[k+c];
            k += 2*blockDim.y;
         }

         k = size-1-2*threadIdx.y;
         for(int i = size/2+threadIdx.y; i < size; i+=blockDim.y)
         {
            out[i+c] = in[k+c];
            k -= 2*blockDim.y;
         }
      }
   }

    __global__ void buildDCTPhases(double* sn, double* cs, const int fftSize, const int size)
    {
       double a = acos(-1.0)/static_cast<double>(fftSize);
       for(int i = threadIdx.x; i < size; i+=blockDim.x)
       {
         double phase = a*i;
         sn[i] = sin(phase);
         cs[i] = cos(phase);
       }
    }

    __global__ void buildDCT4Phases(double * sn, double * cs, const int fftSize, const int size)
    {
       double a = acos(-1.0)/static_cast<double>(fftSize);
       for(int i = threadIdx.x; i < size; i+=blockDim.x)
       {
         double phase = a*(2.0*i+0.5);
         cs[i] = cos(phase);
         sn[i] = sin(phase);
       }
    }

    __global__ void buildDCTOutput(cufftDoubleReal* out, const cufftDoubleComplex *in, const int size, const int blockSize)
    {
       double a = acos(-1.0)/(2.0*size);
       int inSize = size/2+1;
       int imax = (size+1)/2;
       for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
       {
          int ci = j*inSize;
          int co = j*size;
          for(int i = threadIdx.y; i < imax; i+=blockDim.y)
          {
             cufftDoubleReal phase = a*i;
             out[i+co] = cos(phase)*in[i+ci].x + sin(phase)*in[i+ci].y;
          }
          int k = inSize-1-threadIdx.y;
          for(int i = imax+threadIdx.y; i < size; i+=blockDim.y,k-=blockDim.y)
          {
             cufftDoubleReal phase = a*i;
             out[i+co] = cos(phase)*in[k+ci].x - sin(phase)*in[k+ci].y;
          }
       }
    }

    __global__ void buildDCTOutput(cufftDoubleReal* out, const cufftDoubleComplex *in, const int size, const int blockSize, const double* sn, const double* cs)
    {
       int inSize = size/2+1;
       int imax = (size+1)/2;
       for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
       {
          int ci = j*inSize;
          int co = j*size;
          for(int i = threadIdx.y; i < imax; i+=blockDim.y)
          {
             out[i+co] = cs[i]*in[i+ci].x + sn[i]*in[i+ci].y;
          }
          int k = inSize-1-threadIdx.y;
          for(int i = imax+threadIdx.y; i < size; i+=blockDim.y,k-=blockDim.y)
          {
             out[i+co] = cs[i]*in[k+ci].x - sn[i]*in[k+ci].y;
          }
       }
    }

    __global__ void buildIDCTInput(cufftDoubleComplex* out, const cufftDoubleReal *in, const int size, const int blockSize)
    {
       double a = acos(-1.0)/(2.0*size);
       int outSize = size/2+1;
       for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
       {
          int co = j*outSize;
          int ci = j*size;
          int k = size-threadIdx.y;
          int i0 = threadIdx.y;
          if(i0 == 0)
          {
             out[co].x = in[ci];
             out[co].y = 0.0;
             i0 += blockDim.y;
             k -= blockDim.y;
          }
          for(int i = i0; i < outSize; i+=blockDim.y,k-=blockDim.y)
          {
             cufftDoubleReal phase = a*i;
             double cs = cos(phase);
             double sn = sin(phase);
             out[i+co].x = cs*in[i+ci] + sn*in[k+ci];
             out[i+co].y = sn*in[i+ci] - cs*in[k+ci];
          }
       }
    }

    __global__ void buildIDCTInput(cufftDoubleComplex* out, const cufftDoubleReal *in, const int size, const int blockSize, const double *sn, const double *cs)
    {
       int outSize = size/2+1;
       for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
       {
          int co = j*outSize;
          int ci = j*size;
          int k = size-threadIdx.y;
          int i0 = threadIdx.y;
          if(i0 == 0)
          {
             out[co].x = in[ci];
             out[co].y = 0.0;
             i0 += blockDim.y;
             k -= blockDim.y;
          }
          for(int i = i0; i < outSize; i+=blockDim.y,k-=blockDim.y)
          {
             out[i+co].x = cs[i]*in[i+ci] + sn[i]*in[k+ci];
             out[i+co].y = sn[i]*in[i+ci] - cs[i]*in[k+ci];
          }
       }
    }

    __global__ void buildIDCTOutput(cufftDoubleReal* out, const cufftDoubleReal *in, const int size, const int blockSize)
    {
       for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
       {
          int c = j*size;
          int k = 2*threadIdx.y;
          for(int i = threadIdx.y; i < size/2; i+=blockDim.y,k+=2*blockDim.y)
          {
             out[k+c] = in[i+c];
          }

          k = size-1-2*threadIdx.y;
          for(int i = size/2+threadIdx.y; i < size; i+=blockDim.y,k-=2*blockDim.y)
          {
             out[k+c] = in[i+c];
          }
       }
    }

    __global__ void buildDCT4Input(cufftDoubleComplex* out, const cufftDoubleReal *in, const int size, const int blockSize)
    {
       double a = acos(-1.0)/size;
       int outSize = size/2;
       for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
       {
          int co = j*outSize;
          int ci = j*size;
          int k = size-1-2*threadIdx.y;
          for(int i = threadIdx.y; i < outSize; i+=blockDim.y, k-=2*blockDim.y)
          {
             cufftDoubleReal phase = a*i;
             double cs = cos(phase);
             double sn = sin(phase);
             out[i+co].x = cs*in[2*i+ci] + sn*in[k+ci];
             out[i+co].y = -sn*in[2*i+ci] + cs*in[k+ci];
          }
       }
    }

    __global__ void buildDCT4Input(cufftDoubleComplex* out, const cufftDoubleReal *in, const int size, const int blockSize, const double* sn, const double *cs)
    {
       int outSize = size/2;
       for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
       {
          int co = j*outSize;
          int ci = j*size;
          int k = size-1-2*threadIdx.y;
          for(int i = threadIdx.y; i < outSize; i+=blockDim.y, k-=2*blockDim.y)
          {
             out[i+co].x = cs[i]*in[2*i+ci] + sn[i]*in[k+ci];
             out[i+co].y = -sn[i]*in[2*i+ci] + cs[i]*in[k+ci];
          }
       }
    }

    __global__ void buildDCT4Output(cufftDoubleReal* out, const cufftDoubleComplex *in, const int size, const int blockSize, const double scale)
    {
       double a = acos(-1.0)/(2.0*size);
       int inSize = size/2;
       for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
       {
          int co = j*size;
          int ci = j*inSize;
          int k = size/2-1-threadIdx.y;
          for(int i = threadIdx.y; i < inSize; i+=blockDim.y, k-=blockDim.y)
          {
             cufftDoubleReal phase = a*(2.0*i+0.5);
             double cs = cos(phase);
             double sn = sin(phase);
             out[2*i+co] = scale*(cs*in[i+ci].x + sn*in[i+ci].y);
             phase = a*(2.0*k+0.5);
             cs = cos(phase);
             sn = sin(phase);
             out[2*i+1+co] = scale*(sn*in[k+ci].x - cs*in[k+ci].y);
          }
       }
    }

    __global__ void buildDCT4Output(cufftDoubleReal* out, const cufftDoubleComplex *in, const int size, const int blockSize, const double scale, const double *sn, const double *cs)
    {
       int inSize = size/2;
       for(int j = threadIdx.x; j < blockSize; j+=blockDim.x)
       {
          int co = j*size;
          int ci = j*inSize;
          int k = size/2-1-threadIdx.y;
          for(int i = threadIdx.y; i < inSize; i+=blockDim.y, k-=blockDim.y)
          {
             out[2*i+co] = scale*(cs[i]*in[i+ci].x + sn[i]*in[i+ci].y);
             out[2*i+1+co] = scale*(sn[k]*in[k+ci].x - cs[k]*in[k+ci].y);
          }
       }
    }

}
}
}
}
}
}
