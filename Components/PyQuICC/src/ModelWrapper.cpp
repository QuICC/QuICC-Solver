/** 
 * @file ModelWrapper.cpp
 * @brief Source of the model Python interpreter wrapper
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
#include "QuICC/PyQuICC/ModelWrapper.hpp"

// Project includes
//
#include "QuICC/PyQuICC/Config.hpp"
#include "QuICC/PyQuICC/Coordinator.hpp"

namespace QuICC {

namespace PyQuICC {

   ModelWrapper::ModelWrapper()
      : mpModule(NULL), mpFunc(NULL), mpModel(NULL), mpMethod(NULL)
   {
      this->init();
   }

   ModelWrapper::~ModelWrapper()
   {
      this->finalize();
   }

   void ModelWrapper::init()
   {
      Coordinator::init();
   }

   void ModelWrapper::import(const std::string& module)
   {
      // Cleanup
      if(this->mpModule != NULL)
      {
         Py_CLEAR(this->mpModule);
      }

      // Get string object for module name
      PyObject* pName;
      pName = PyUnicode_FromString(("quicc.model."+module).c_str());

      // Import module
      this->mpModule = PyImport_Import(pName);

      // Release pName
      Py_DECREF(pName);

      if(this->mpModule == NULL)
      {
         PyErr_Print();
         throw std::logic_error("Python module import error!");
      }
   }

   void ModelWrapper::createModel(const std::string& model)
   {
      // Cleanup
      if(this->mpModel != NULL)
      {
         Py_CLEAR(this->mpModel);
      }

      // Create model object

      this->setFunction(model);
      this->mpModel = this->callFunction();

      if(this->mpModel == NULL)
      {
         PyErr_Print();
         throw std::logic_error("Python model creation error!");
      }
   }

   void ModelWrapper::enableGalerkin(const bool flag)
   {
      // Set the galerkin flag for the python dispatchers
      if(flag)
      {
         PyObject_SetAttrString(this->mpModel, (char *)"use_galerkin", Py_True);
      } else
      {
         PyObject_SetAttrString(this->mpModel, (char *)"use_galerkin", Py_False);
      }
   }

   void ModelWrapper::createModel(const std::string& model, const std::string& specialization)
   {
      // Check module is setup
      this->checkModule();

      bool hasSpecialization = PyObject_HasAttrString(this->mpModule, (model+specialization).c_str());
      if(hasSpecialization)
      {
         this->createModel(model + specialization);
      } else
      {
         this->createModel(model);
      }
   }

   void ModelWrapper::setFunction(const std::string& func)
   {
      // Check module is setup
      this->checkModule();

      // Cleanup
      if(this->mpFunc != NULL)
      {
         Py_CLEAR(this->mpFunc);
      }

      // Get Python function object
      this->mpFunc = PyObject_GetAttrString(this->mpModule, func.c_str());

      // Check for successfully loading function
      if(! (this->mpFunc && PyCallable_Check(this->mpFunc)))
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
         throw std::logic_error("Python function loading error!");
      }
   }

   void ModelWrapper::setFunction(const std::string& func, const std::string& submod)
   {
      // Check module is setup
      this->checkModule();

      // Cleanup
      if(this->mpFunc != NULL)
      {
         Py_CLEAR(this->mpFunc);
      }

      // Get Python function object
      PyObject *pTmp = PyObject_GetAttrString(this->mpModule, submod.c_str());
      if(pTmp)
      {
         this->mpFunc = PyObject_GetAttrString(pTmp, func.c_str());
         Py_DECREF(pTmp);

         // Check for successfully loading function
         if(! (this->mpFunc && PyCallable_Check(this->mpFunc)))
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

   void ModelWrapper::setMethod(const std::string& method)
   {
      // Check model is setup
      this->checkModel();

      // Cleanup
      if(this->mpMethod != NULL)
      {
         Py_CLEAR(this->mpMethod);
      }

      // Get model method object
      this->mpMethod = PyObject_GetAttrString(this->mpModel, method.c_str());

      // Check for successfully loaded method
      if(! (this->mpMethod && PyCallable_Check(this->mpMethod)))
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
         throw std::logic_error("Python model method loading error!");
      }
   }

   PyObject* ModelWrapper::callFunction()
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(this->mpFunc, NULL);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw std::logic_error("Python function call error!");
      }

      return pValue;
   }

   PyObject* ModelWrapper::callFunction(PyObject* pArgs)
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(this->mpFunc, pArgs);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw std::logic_error("Python function with arguments call error!");
      }

      return pValue;
   }

   PyObject* ModelWrapper::callMethod()
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(this->mpMethod, NULL);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw std::logic_error("Python method call error!");
      }

      return pValue;
   }

   PyObject* ModelWrapper::callMethod(PyObject* pArgs)
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(this->mpMethod, pArgs);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw std::logic_error("Python method with arguments call error!");
      }

      return pValue;
   }

   void ModelWrapper::fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat) 
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
      this->setFunction("triplets", "utils");
      pValue = this->callFunction(pArgs);
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

   void ModelWrapper::fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat) 
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
      this->setFunction("triplets", "utils");
      pValue = this->callFunction(pArgs);
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

   void ModelWrapper::fillMatrix(Matrix& rMatrix, PyObject* pPyMat)
   {
      PyObject *pArgs, *pValue;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));
      Py_DECREF(pValue);

      // Allocate matrix
      rMatrix.resize(rows,cols);

      // Convert Python matrix into a list
      pArgs = PyTuple_New(1);
      Py_INCREF(pPyMat);
      PyTuple_SetItem(pArgs, 0, pPyMat);
      this->setFunction("nparr2list", "utils");
      pValue = this->callFunction(pArgs);
      Py_DECREF(pArgs);

      long int count = 0;
      for(int j = 0; j < cols; j++)
      {
         for (int i = 0; i < rows; i++)
         {
            rMatrix(i, j) = PyFloat_AsDouble(PyList_GetItem(pValue, count));
            count += 1;
         }
      }
      Py_DECREF(pValue);
   }

   void ModelWrapper::fillArray(Array& rArray, PyObject* pPyArray)
   {
      PyObject *pArgs, *pValue;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyArray, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      Py_DECREF(pValue);

      // Allocate matrix
      rArray.resize(rows);

      // Convert Python matrix into a list
      pArgs = PyTuple_New(1);
      Py_INCREF(pPyArray);
      PyTuple_SetItem(pArgs, 0, pPyArray);
      this->setFunction("nparr2list", "utils");
      pValue = this->callFunction(pArgs);
      Py_DECREF(pArgs);

      long int count = 0;
      for (int i = 0; i < rows; i++)
      {
         rArray(i) = PyFloat_AsDouble(PyList_GetItem(pValue, count));
         count += 1;
      }

      Py_DECREF(pValue);
   }

   void ModelWrapper::checkModule()
   {
      if(this->mpModule == NULL)
      {
         throw std::logic_error("Python module is not setup!");
      }
   }

   void ModelWrapper::checkModel()
   {
      // Check model is setup
      if(this->mpModel == NULL)
      {
         throw std::logic_error("Python model is not setup!");
      }
   }

   void ModelWrapper::cleanup()
   {
      // Clean up
      if(this->mpFunc != NULL)
      {
         Py_CLEAR(this->mpFunc);
      }
      if(this->mpMethod != NULL)
      {
         Py_CLEAR(this->mpMethod);
      }
   }

   void ModelWrapper::finalize()
   {
      // Clean up
      this->cleanup();

      // Clear model and module
      if(this->mpModel != NULL)
      {
         Py_CLEAR(this->mpModel);
      }
      if(this->mpModule != NULL)
      {
         Py_CLEAR(this->mpModule);
      }

      // Finalize
      Coordinator::finalize();
   }

}
}
