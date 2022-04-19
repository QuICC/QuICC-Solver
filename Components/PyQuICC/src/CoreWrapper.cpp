/**
 * @file CoreWrapper.cpp
 * @brief Source of Python interpreter wrapper
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
#include "QuICC/PyQuICC/CoreWrapper.hpp"

// Project includes
//
#include "QuICC/PyQuICC/Config.hpp"
#include "QuICC/PyQuICC/Coordinator.hpp"

namespace QuICC {

namespace PyQuICC {

   PyObject* CoreWrapper::mpModule = NULL;

   PyObject* CoreWrapper::mpFunc = NULL;

   PyObject* CoreWrapper::mpClass = NULL;

   PyObject* CoreWrapper::mpMethod = NULL;

   void CoreWrapper::init()
   {
      Coordinator::init();
   }

   void CoreWrapper::import(const std::string& module)
   {
      // Cleanup
      if(CoreWrapper::mpModule != NULL)
      {
         Py_CLEAR(CoreWrapper::mpModule);
      }

      // Get string object for module name
      PyObject* pName;
      pName = PyUnicode_FromString(module.c_str());

      // Import module
      CoreWrapper::mpModule = PyImport_Import(pName);

      // Release pName
      Py_DECREF(pName);

      if(CoreWrapper::mpModule == NULL)
      {
         PyErr_Print();
         throw std::logic_error("Python module import error!");
      }
   }

   void CoreWrapper::createClass(const std::string& name)
   {
      // Cleanup
      if(CoreWrapper::mpClass != NULL)
      {
         Py_CLEAR(CoreWrapper::mpClass);
      }

      // Create class object

      CoreWrapper::setFunction(name);
      CoreWrapper::mpClass = CoreWrapper::callFunction();

      if(CoreWrapper::mpClass == NULL)
      {
         PyErr_Print();
         throw std::logic_error("Python class object creation error!");
      }
   }

   void CoreWrapper::setFunction(const std::string& func)
   {
      // Cleanup
      if(CoreWrapper::mpFunc != NULL)
      {
         Py_CLEAR(CoreWrapper::mpFunc);
      }

      // Get Python function object
      CoreWrapper::mpFunc = PyObject_GetAttrString(CoreWrapper::mpModule, func.c_str());

      // Check for successfully loading function
      if(! (CoreWrapper::mpFunc && PyCallable_Check(CoreWrapper::mpFunc)))
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
         throw std::logic_error("Python function loading error!");
      }
   }

   void CoreWrapper::setFunction(const std::string& func, const std::string& submod)
   {
      // Cleanup
      if(CoreWrapper::mpFunc != NULL)
      {
         Py_CLEAR(CoreWrapper::mpFunc);
      }

      // Get Python function object
      PyObject *pTmp = PyObject_GetAttrString(CoreWrapper::mpModule, submod.c_str());

      if(pTmp)
      {
         CoreWrapper::mpFunc = PyObject_GetAttrString(pTmp, func.c_str());
         Py_DECREF(pTmp);

         // Check for successfully loading function
         if(! (CoreWrapper::mpFunc && PyCallable_Check(CoreWrapper::mpFunc)))
         {
            if(PyErr_Occurred())
            {
               PyErr_Print();
            }
            throw std::logic_error("Python function loading error!");
         }
      } else
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
      }
   }

   void CoreWrapper::setMethod(const std::string& method)
   {
      // Cleanup
      if(CoreWrapper::mpMethod != NULL)
      {
         Py_CLEAR(CoreWrapper::mpMethod);
      }

      // Get class method object
      CoreWrapper::mpMethod = PyObject_GetAttrString(CoreWrapper::mpClass, method.c_str());

      // Check for successfully loaded method
      if(! (CoreWrapper::mpMethod && PyCallable_Check(CoreWrapper::mpMethod)))
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
         throw std::logic_error("Python class method loading error!");
      }
   }

   PyObject* CoreWrapper::callFunction()
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(CoreWrapper::mpFunc, NULL);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw std::logic_error("Python function call error!");
      }

      return pValue;
   }

   PyObject* CoreWrapper::callFunction(PyObject* pArgs)
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(CoreWrapper::mpFunc, pArgs);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw std::logic_error("Python function with arguments call error!");
      }

      return pValue;
   }

   PyObject* CoreWrapper::callMethod()
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(CoreWrapper::mpMethod, NULL);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw std::logic_error("Python method call error!");
      }

      return pValue;
   }

   PyObject* CoreWrapper::callMethod(PyObject* pArgs)
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(CoreWrapper::mpMethod, pArgs);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw std::logic_error("Python method with arguments call error!");
      }

      return pValue;
   }

   void CoreWrapper::fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat)
   {
      PyObject *pArgs, *pValue, *pTmp;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));
      Py_DECREF(pValue);

      // Convert Python matrix into triplets
      pArgs = PyTuple_New(1);
      Py_INCREF(pPyMat);
      PyTuple_SetItem(pArgs, 0, pPyMat);
      CoreWrapper::setFunction("triplets", "utils");
      pValue = CoreWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      long int len = PyList_Size(pValue);
      std::vector<Triplet> triplets;
      triplets.reserve(len);
      long int row;
      long int col;
      MHDFloat val;
      for(int i = 0; i < len; i++)
      {
         pTmp = PyList_GetItem(pValue, i);
         row = PyLong_AsLong(PyTuple_GetItem(pTmp,0));
         col = PyLong_AsLong(PyTuple_GetItem(pTmp,1));
         val = PyFloat_AsDouble(PyTuple_GetItem(pTmp,2));
         triplets.push_back(Triplet(row,col,val));
      }
      Py_DECREF(pValue);

      // Build matrix
      rMatrix.resize(rows,cols);
      rMatrix.setFromTriplets(triplets.begin(), triplets.end());
   }

   void CoreWrapper::fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat)
   {
      PyObject *pArgs, *pValue, *pTmp;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));
      Py_DECREF(pValue);

      // Convert Python matrix into triplets
      pArgs = PyTuple_New(1);
      Py_INCREF(pPyMat);
      PyTuple_SetItem(pArgs, 0, pPyMat);
      CoreWrapper::setFunction("triplets", "utils");
      pValue = CoreWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      long int len = PyList_Size(pValue);
      std::vector<Triplet> realTriplets;
      std::vector<Triplet> imagTriplets;
      realTriplets.reserve(len);
      imagTriplets.reserve(len);
      long int row;
      long int col;
      MHDFloat val;
      for(int i = 0; i < len; i++)
      {
         pTmp = PyList_GetItem(pValue, i);
         row = PyLong_AsLong(PyTuple_GetItem(pTmp,0));
         col = PyLong_AsLong(PyTuple_GetItem(pTmp,1));
         if(PyComplex_Check(PyTuple_GetItem(pTmp,2)))
         {
            val = PyComplex_RealAsDouble(PyTuple_GetItem(pTmp,2));
            realTriplets.push_back(Triplet(row,col,val));
            val = PyComplex_ImagAsDouble(PyTuple_GetItem(pTmp,2));
            imagTriplets.push_back(Triplet(row,col,val));
         } else
         {
            val = PyFloat_AsDouble(PyTuple_GetItem(pTmp,2));
            realTriplets.push_back(Triplet(row,col,val));
         }
      }
      Py_DECREF(pValue);

      // Build matrix
      rMatrix.real().resize(rows,cols);
      rMatrix.imag().resize(rows,cols);
      rMatrix.real().setFromTriplets(realTriplets.begin(), realTriplets.end());

      if(imagTriplets.size() > 0)
      {
         rMatrix.imag().setFromTriplets(imagTriplets.begin(), imagTriplets.end());
      }
   }

   void CoreWrapper::cleanup()
   {
      // Clean up
      if(CoreWrapper::mpFunc != NULL)
      {
         Py_CLEAR(CoreWrapper::mpFunc);
      }
      if(CoreWrapper::mpMethod != NULL)
      {
         Py_CLEAR(CoreWrapper::mpMethod);
      }
   }

   void CoreWrapper::finalize()
   {
      // Clean up
      CoreWrapper::cleanup();

      // Clear class and module
      if(CoreWrapper::mpClass != NULL)
      {
         Py_CLEAR(CoreWrapper::mpClass);
      }
      if(CoreWrapper::mpModule != NULL)
      {
         Py_CLEAR(CoreWrapper::mpModule);
      }

      // Finalize
      Coordinator::finalize();
   }

   CoreWrapper::CoreWrapper()
   {
   }

   CoreWrapper::~CoreWrapper()
   {
   }

}
}
