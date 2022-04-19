/**
 * @file Pardiso_Real.cpp
 * @brief Implementation of the C++ wrapper functions for the real valued Pardiso routines 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#include "../External/Interfaces/Pardiso_Interface.hpp"

extern "C"
{
   /**
    * @brief Declaration for the real valued pardiso routine
    */
   void pardiso(void   *, int *,   int *, int *, int *, int *, 
                double *, int *,   int *, int *,   int *, int *,
                int *, double *, double *, int *, double *);
}

namespace Eigen {

namespace Pardiso {

   void pardisoinit(void *pt, int *solver, int *iparm, double *dparm, int *error, double)
   {
      int mtype = 11;
      ::pardisoinit(pt, &mtype, solver, iparm, dparm, error);
   } 

   void pardiso(void *pt, const int *maxfct, const int *mnum, const int *phase, const int *n, 
               const double *a, const int *ia, const int *ja, const int *perm, const int *nrhs, int *iparm,
               const int *msglvl, double *b, double *x, int *error, double *dparm, double)
   {
      int mtype = 11;
      ::pardiso(pt, const_cast<int *>(maxfct), const_cast<int *>(mnum), &mtype, const_cast<int *>(phase), const_cast<int *>(n), const_cast<double *>(a), const_cast<int *>(ia), const_cast<int *>(ja), const_cast<int *>(perm), const_cast<int *>(nrhs), iparm, const_cast<int *>(msglvl), b, x, error, dparm);
   }

} // end namespace Pardiso

} // end namespace Eigen
