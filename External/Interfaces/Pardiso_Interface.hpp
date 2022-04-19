/**
 * @file Pardiso_Interface.hpp
 * @brief Declarations needed to use the Pardiso routines in the C++ 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PARDISO_INTERFACE_HPP
#define PARDISO_INTERFACE_HPP

extern "C"
{
   /**
    * @brief Declaration for the pardisoinit routine
    */
   void pardisoinit(void *, int *, int *, int *, double *, int *);

   /**
    * @brief Declaration for the real valued pardiso_chkmatrix routine
    */
   void pardiso_chkmatrix(int *, int *, double *, int *, int *, int *);

   /**
    * @brief Declaration for the real valued pardiso_chkvec routine
    */
   void pardiso_chkvec(int *, int *, double *, int *);

   /**
    * @brief Declaration for the real valued pardiso_printstats routine
    */
   void pardiso_printstats(int *, int *, double *, int *, int *, int *, double *, int *);

   /**
    * @brief Typedef for the complex values struct used by pardiso
    */
   typedef struct { double re; double i;} pardiso_complex;

   /**
    * @brief Declaration for the complex valued pardiso_chkmatrix routine
    */
   void pardiso_chkmatrix_z(int *, int *, pardiso_complex *, int *, int *, int *);

   /**
    * @brief Declaration for the complex valued pardiso_chkvec routine
    */
   void pardiso_chkvec_z(int *, int *, pardiso_complex *, int *);

   /**
    * @brief Declaration for the complex valued pardiso_printstats routine
    */
   void pardiso_printstats_z(int *, int *, pardiso_complex *, int *, int *, int *, pardiso_complex *, int *);
}

#endif // PARDISO_INTERFACE_HPP
