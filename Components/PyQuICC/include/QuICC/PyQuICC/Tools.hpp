/**
 * @file Tools.hpp
 * @brief Tools to work with Python objects
 */

#ifndef QUICC_PYQUICC_TOOLS_HPP
#define QUICC_PYQUICC_TOOLS_HPP

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
#include "QuICC/NonDimensional/INumber.hpp"

namespace QuICC {

namespace PyQuICC {

   /**
    * @brief Tools to work with Python objects
    */
   class Tools
   {
      public:
         /**
          * @brief Make a tuple from integer array
          */
         static PyObject* makeTuple(const ArrayI& val);

         /**
          * @brief Make a tuple from vector of double
          */
         static PyObject* makeTuple(const std::vector<MHDFloat>& val);

         /**
          * @brief Make a list from vector of int
          */
         static PyObject* makeList(const std::vector<int>& val);

         /**
          * @brief Make a list of list from vector of vector of int
          */
         static PyObject* makeList(const std::vector<std::vector<int> >& val);

         /**
          * @brief Make a dictionary
          */
         template <typename T> static PyObject* makeDict(const T& map);

         /**
          * @brief Make a dictionary
          */
         template <typename T1, typename T2> static PyObject* makeDict(const T1& key, const T2& val);

         /**
          * @brief Append to dictionary
          */
         static void addDict(const std::vector<std::string>& key, const std::vector<MHDFloat>& val, PyObject* pDict);

         /**
          * @brief Append to dictionary
          */
         static void addDict(const std::map<std::string,int>& map, PyObject* pDict);

         /**
          * @brief Append to dictionary
          */
         static void addDict(const std::map<std::size_t,NonDimensional::SharedINumber>& map, PyObject* pDict);

         /**
          * @brief Append to dictionary
          */
         static void addDict(const std::map<std::string,MHDFloat>& map, PyObject* pDict);

         /**
          * @brief Append to dictionary
          */
         static void addDict(const std::map<std::string,std::string>& map, PyObject *pDict);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<bool> &rList, PyObject *pList);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<std::size_t> &rList, PyObject *pList);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<std::string> &rList, PyObject *pList);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<std::pair<std::size_t,FieldComponents::Spectral::Id> > &rList, PyObject *pList);

         /**
          * @brief Get data from dict
          */
         static void getDict(std::map<std::size_t,NonDimensional::SharedINumber> &rMap, PyObject *pDict, const bool replace);

         /**
          * @brief Get data from dict
          */
         static void getDict(std::map<std::string,MHDFloat> &rMap, PyObject *pDict, const bool replace);

         /**
          * @brief Convert NumPy array to list
          */
         static PyObject* nparr2List(PyObject *pArr);

         /**
          * @brief Convert Scipy sparse matrix to triplets
          */
         static PyObject* sparse2triplets(PyObject *pSpMat);

         /**
          * @brief Fill Array with NumPy array data
          */
         static void fillArray(Array& rArray, PyObject* pPyArr);

         /**
          * @brief Fill Matrix with NumPy 2D array data
          */
         static void fillMatrix(Matrix& rMatrix, PyObject* pPyMat);

         /**
          * @brief Fill sparse matrix with data from Python call
          */
         static void fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat);

         /**
          * @brief Fill sparse matrix with data from Python call
          */
         static void fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat);
  
         /**
          * @brief Cleanup wrapper without finalize
          */
         static void cleanup();
  
         /**
          * @brief Finalise the Python interpreter
          */
         static void finalize();
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         Tools() = default;

         /**
          * @brief Destructor
          */
         ~Tools() = default;
   };

   template <typename T> PyObject* Tools::makeDict(const T& map)
   {
      PyObject *pDict;
      pDict = PyDict_New();

      Tools::addDict(map, pDict);

      return pDict;
   }

   template <typename T1, typename T2> PyObject* Tools::makeDict(const T1& key, const T2& val)
   {
      PyObject *pDict;
      pDict = PyDict_New();

      Tools::addDict(key, val, pDict);

      return pDict;
   }

}
}

#endif // QUICC_PYQUICC_TOOLS_HPP
