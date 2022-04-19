/** 
 * @file NumpyWrapper.cpp
 * @brief Source of Python interpreter wrapper for Numpy types
 */

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/PyQuICC/NumpyWrapper.hpp"

// Project includes
//

namespace QuICC {

namespace PyQuICC {

   PyObject* NumpyWrapper::mpModule = NULL;

   PyObject* NumpyWrapper::mpFunc = NULL;

   PyObject* NumpyWrapper::mpClass = NULL;

   PyObject* NumpyWrapper::mpMethod = NULL;

   void NumpyWrapper::init()
   {
      if(_import_array() < 0)
      {
         throw std::logic_error("Initialization of NumPy C-API failed");
      }

      CoreWrapper::init();
   }

   void NumpyWrapper::getMatrix(Matrix& rMatrix, PyArrayObject* pMat)
   {
         // TODO: some precondition-checking on the number of dimensions and the size
         // get the size of Python matrix
         long* dims;
         dims = PyArray_DIMS(pMat);
         //dims = ((PyArrayObject_fields *)pMat)->dimensions;

         // resize the matrix to correct size and assign the pointer
         // note that that numpy default storage is RowMajor whereas Eigen 3.3.1 is ColumnMajor
         //rMatrix = Eigen::Map<Eigen::Matrix<MHDFloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >((double*)PyArray_DATA(pMat),dims[0],dims[1]);
         rMatrix = Eigen::Map<Matrix>((double*)PyArray_DATA(pMat),dims[0],dims[1]);
         //rMatrix = Eigen::Map<Matrix>( ( double*) ((PyArrayObject_fields *)pMat)->data, dims[0], dims[1]);

   }

   void NumpyWrapper::getVector(Array& rVector, PyArrayObject* pVec)
   {

         // get the  lenght of the vector
         int len = PyArray_DIM(pVec,0);

         // resize the eigen::vector and assign the pointer
         rVector = Eigen::Map<Eigen::VectorXd>((double*)PyArray_DATA(pVec),len);
   }

   void NumpyWrapper::getVector(ArrayZ& rVector, PyArrayObject* pVec)
   {

         // get the  lenght of the vector
         int len = PyArray_DIM(pVec,0);

         // resize the eigen::vector and assign the pointer
         rVector = Eigen::Map<Eigen::VectorXcd>((std::complex<double>*)PyArray_DATA(pVec),len);
   }

   NumpyWrapper::NumpyWrapper()
      :CoreWrapper()
   {

   }

   NumpyWrapper::~NumpyWrapper()
   {
   }

}
}
