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
#include "QuICC/PyQuICC/Tools.hpp"

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

   void CoreWrapper::fillArray(Array& rArray, PyObject* pPyArray)
   {
      Tools::fillArray(rArray, pPyArray);
   }

   void CoreWrapper::fillMatrix(Matrix& rMatrix, PyObject* pPyMat)
   {
      Tools::fillMatrix(rMatrix, pPyMat);
   }

   void CoreWrapper::fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat)
   {
      Tools::fillMatrix(rMatrix, pPyMat);
   }

   void CoreWrapper::fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat)
   {
      Tools::fillMatrix(rMatrix, pPyMat);
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
