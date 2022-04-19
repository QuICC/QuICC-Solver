/**
 * @file Pardiso_Complex.cpp
 * @brief Declarations needed to use the complex Pardiso routines in the C++ 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "../External/Interfaces/Pardiso_Interface.hpp"
#include <complex>

extern "C"
{
   /**
    * @brief Declaration for the complex valued pardiso routine
    */
   void pardiso(void *, int *, int *, int *, int *, int *, 
                pardiso_complex *, int *,   int *, int *,   int *, int *,
                int *, pardiso_complex *, pardiso_complex *, int *, double *);
}

namespace Eigen {

namespace Pardiso {

   void pardisoinit(void *pt, int *solver, int *iparm, double *dparm, int *error, std::complex<double>)
   {
      int mtype = 13;
      ::pardisoinit(pt, &mtype, solver, iparm, dparm, error);
   }

   void pardiso(void *pt, const int *maxfct, const int *mnum, const int *phase, const int *n, 
               const std::complex<double> *a, const int *ia, const int *ja, const int *perm, const int *nrhs, int *iparm,
               const int *msglvl, std::complex<double> *b, std::complex<double> *x, int *error, double *dparm, std::complex<double>)
   {
      int mtype = 13;
      ::pardiso(pt, const_cast<int *>(maxfct), const_cast<int *>(mnum), &mtype, const_cast<int *>(phase), const_cast<int *>(n), const_cast<pardiso_complex *>(reinterpret_cast<const pardiso_complex *>(a)), const_cast<int *>(ia), const_cast<int *>(ja), const_cast<int *>(perm), const_cast<int *>(nrhs), iparm, const_cast<int *>(msglvl), reinterpret_cast<pardiso_complex *>(b), reinterpret_cast<pardiso_complex *>(x), error, dparm);
   }

} // end namespace Pardiso

} // end namespace Eigen
