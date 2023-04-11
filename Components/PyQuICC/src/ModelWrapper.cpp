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
#include "QuICC/PyQuICC/Tools.hpp"

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
      pName = PyUnicode_FromString(("quicc_solver.model."+module).c_str());

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

   void ModelWrapper::enableSplitEquation(const bool flag)
   {
      // Split equation are not implemented in Python backend
      if(flag)
      {
         throw std::logic_error("Python backend does not implemented the split equations");
      }

//      // Set the split equation flag for the python dispatchers
//      if(flag)
//      {
//         PyObject_SetAttrString(this->mpModel, (char *)"use_splitequation", Py_True);
//      } else
//      {
//         PyObject_SetAttrString(this->mpModel, (char *)"use_splitequation", Py_False);
//      }
   }

   void ModelWrapper::enableLinearized(const bool flag)
   {
      // Set the galerkin flag for the python dispatchers
      if(flag)
      {
         PyObject_SetAttrString(this->mpModel, (char *)"linearize", Py_True);
      } else
      {
         PyObject_SetAttrString(this->mpModel, (char *)"linearize", Py_False);
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

   void ModelWrapper::fillArray(Array& rArray, PyObject* pPyArray)
   {
      Tools::fillArray(rArray, pPyArray);
   }

   void ModelWrapper::fillMatrix(Matrix& rMatrix, PyObject* pPyMat)
   {
      Tools::fillMatrix(rMatrix, pPyMat);
   }

   void ModelWrapper::fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat) 
   {
      Tools::fillMatrix(rMatrix, pPyMat);
   }

   void ModelWrapper::fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat) 
   {
      Tools::fillMatrix(rMatrix, pPyMat);
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
