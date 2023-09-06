/**
 * @file ModelWrapper.hpp
 * @brief Static class wrapper to the Python model embedding
 */

#ifndef QUICC_PYQUICC_MODELWRAPPER_HPP
#define QUICC_PYQUICC_MODELWRAPPER_HPP

// First include
//
#include "QuICC/PyQuICC/SystemHeader.hpp"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace PyQuICC {

   /**
    * @brief Wrapper to the Python model embedding
    */
   class ModelWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         ModelWrapper();

         /**
          * @brief Destructor
          */
         ~ModelWrapper();

         /**
          * @brief Enable Galerkin?
          */
         void enableGalerkin(const bool flag);

         /**
          * @brief Enable split equation?
          */
         void enableSplitEquation(const bool flag);

         /**
          * @brief Enable linearized equation?
          */
         void enableLinearized(const bool flag);

         /**
          * @brief Import the Python module
          */
         void import(const std::string& module);

         /**
          * @brief Create the Python model class
          */
         void createModel(const std::string& model);

         /**
          * @brief Create the Python model class, first trying model+specialization
          */
         void createModel(const std::string& model, const std::string& specialization);

         /**
          * @brief Set the Python function object
          */
         void setFunction(const std::string& func);

         /**
          * @brief Set the Python function object
          */
         void setFunction(const std::string& func, const std::string& submod);

         /**
          * @brief Call the Python function object
          */
         PyObject* callFunction();

         /**
          * @brief Call the Python function object with arguments
          */
         PyObject* callFunction(PyObject* pArgs);

         /**
          * @brief Set the Python model method object
          */
         void setMethod(const std::string& method);

         /**
          * @brief Call the Python model method object
          */
         PyObject* callMethod();

         /**
          * @brief Call the Python model method object with arguments
          */
         PyObject* callMethod(PyObject* pArgs);

         /**
          * @brief Fill array with data from Python call
          */
         void fillArray(Array& rArray, PyObject* pPyArray);

         /**
          * @brief Fill dense matrix with data from Python call
          */
         void fillMatrix(Matrix& rMatrix, PyObject* pPyMat);

         /**
          * @brief Fill sparse matrix with data from Python call
          */
         void fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat);

         /**
          * @brief Fill sparse matrix with data from Python call
          */
         void fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat);
  
         /**
          * @brief Cleanup wrapper
          */
         void cleanup();
         
      protected:
         /**
          * @brief Initialize the wrapper
          */
         void init();

         /**
          * @brief Finalise the wrapper
          */
         void finalize();

      private:
         /**
          * @brief Check module is not null
          */
         void checkModule();

         /**
          * @brief Check model is not null
          */
         void checkModel();

         /**
          * @brief Python module object
          */
         PyObject* mpModule;

         /**
          * @brief Python function object
          */
         PyObject* mpFunc;

         /**
          * @brief Python model object
          */
         PyObject* mpModel;

         /**
          * @brief Python model method object
          */
         PyObject* mpMethod;
   };

}
}

#endif // QUICC_PYQUICC_MODELWRAPPER_HPP
