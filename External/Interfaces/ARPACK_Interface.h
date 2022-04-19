/**
 * @file ARPACK_Interface.h
 * @brief Declarations needed to use ARPACK routines in the C++
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ARPACK_INTERFACE_HPP
#define ARPACK_INTERFACE_HPP

#include <complex>

extern "C"
{
   /**
    * @brief Declaration for the ZNAUPD routine
    */
   void znaupd_(int* ido, const char* bmat, int* n, const char* which, int* nev, double* tol, std::complex<double>* resid, int* ncv, std::complex<double>* v, int* ldv, int* iparam, int* ipntr, std::complex<double>* workd, std::complex<double>* workl, int* lworkl, double* rwork, int* info);

   /**
    * @brief Declaration for the ZNEUPD routine
    */
   void zneupd_(int* rvec, const char* howmny, int* select, std::complex<double>* d, std::complex<double>* z, int* ldz , std::complex<double>* sigma, std::complex<double>* workev, const char* bmat, int* n, const char* which, int* nev, double* tol, std::complex<double>* resid, int* ncv, std::complex<double>* v, int* ldv, int* iparam, int* ipntr, std::complex<double>* workd, std::complex<double>* workl, int* lworkl, double* rwork, int* info);
}

#endif // ARPACK_INTERFACE_HPP
