/**
 * @file NumpyWrapper.hpp
 * @brief Small wrapper for Python embedding for Numpy arrays
 */

#ifndef QUICC_PYQUICC_NUMPYWRAPPER_HPP
#define QUICC_PYQUICC_NUMPYWRAPPER_HPP


#define NPY_NO_DEPRECATED_API NPY_1_14_API_VERSION

// First include
//
#include "QuICC/PyQuICC/SystemHeader.hpp"

// External include
#include <numpy/ndarrayobject.h>
#include <numpy/npy_3kcompat.h>

// Project includes
//
#include "QuICC/PyQuICC/CoreWrapper.hpp"

namespace QuICC {

namespace PyQuICC {

   class NumpyWrapper: public CoreWrapper
   {

      public:

         static void init();

          /*
           * @brief Fill a full matrix (Eigen::MatrixXd) with data from PyObject
           */
          static void getMatrix(Matrix& rMatrix, PyArrayObject* pMat);

          /*
           * @brief Fill a vector (Eigen::VectorXd) with data from PyObject
           */
          static void getVector(Array& rArray, PyArrayObject* pVec);

          /*
           * @brief Fill a complex vector (Eigen::VectorXcd) with data from PyObjct
           */
          static void getVector(ArrayZ& rArrayZ, PyArrayObject* pVec);

      protected:
          /**
           * @brief Constructor
           */
          NumpyWrapper();

          /**
           * @brief Destructor
           */
          ~NumpyWrapper();

      private:
          /**
           * @brief Python module object
           */
          static PyObject* mpModule;

          /**
           * @brief Python function object
           */
          static PyObject* mpFunc;

          /**
           * @brief Python class object
           */
          static PyObject* mpClass;

          /**
           * @brief Python class method object
           */
          static PyObject* mpMethod;
   };

}
}

#endif // QUICC_PYQUICC_NUMPYWRAPPER_HPP
