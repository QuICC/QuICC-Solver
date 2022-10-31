/** 
 * @file Tools.cpp
 * @brief Source of the tools to work with Python objects
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
#include "QuICC/PyQuICC/Tools.hpp"

// Project includes
//
#include "QuICC/Hasher.hpp"
#include "QuICC/NonDimensional/Coordinator.hpp"
#include "QuICC/Tools/HumanToId.hpp"

namespace QuICC {

namespace PyQuICC {

   PyObject* Tools::makeTuple(const ArrayI& val)
   {
      PyObject *pTuple, *pValue;

      pTuple = PyTuple_New(val.size());
      for(unsigned int i = 0; i < val.size(); i++)
      {
         pValue = PyLong_FromLong(val(i));
         PyTuple_SetItem(pTuple, i, pValue); 
      }

      return pTuple;
   }

   PyObject* Tools::makeTuple(const std::vector<MHDFloat>& val)
   {
      PyObject *pTuple, *pValue;

      pTuple = PyTuple_New(val.size());
      for(unsigned int i = 0; i < val.size(); i++)
      {
         pValue = PyFloat_FromDouble(val.at(i));
         PyTuple_SetItem(pTuple, i, pValue); 
      }

      return pTuple;
   }

   PyObject* Tools::makeList(const std::vector<int>& val)
   {
      PyObject *pList, *pValue;

      pList = PyList_New(val.size());
      for(unsigned int i = 0; i < val.size(); i++)
      {
         pValue = PyLong_FromLong(val.at(i));
         PyList_SetItem(pList, i, pValue); 
      }

      return pList;
   }

   PyObject* Tools::makeList(const std::vector<std::vector<int> >& val)
   {
      PyObject *pList, *pValue;

      pList = PyList_New(val.size());
      for(unsigned int i = 0; i < val.size(); i++)
      {
         pValue = Tools::makeList(val.at(i));
         PyList_SetItem(pList, i, pValue); 
      }

      return pList;
   }

   void Tools::addDict(const std::vector<std::string>& key, const std::vector<MHDFloat>& val, PyObject* pDict)
   {
      if(!PyDict_Check(pDict))
      {
         throw std::logic_error("Python object is not a dictionary. Cannot append");
      }

      PyObject *pKey, *pValue;

      for(unsigned int i = 0; i < key.size(); i++)
      {
         pKey = PyUnicode_FromString(key.at(i).c_str());
         pValue = PyFloat_FromDouble(val.at(i));
         PyDict_SetItem(pDict, pKey, pValue);
         Py_DECREF(pValue);
         Py_DECREF(pKey);
      }
   }

   void Tools::addDict(const std::map<std::string, int>& map, PyObject *pDict)
   {
      if(!PyDict_Check(pDict))
      {
         throw std::logic_error("Python object is not a dictionary. Cannot append");
      }

      PyObject *pKey, *pValue;

      for(auto mapIt = map.cbegin(); mapIt != map.cend(); ++mapIt)
      {
         pKey = PyUnicode_FromString(mapIt->first.c_str());
         pValue = PyLong_FromLong(mapIt->second);
         PyDict_SetItem(pDict, pKey, pValue);
         Py_DECREF(pValue);
         Py_DECREF(pKey);
      }
   }

   void Tools::addDict(const std::map<std::size_t, NonDimensional::SharedINumber>& map, PyObject *pDict)
   {
      if(!PyDict_Check(pDict))
      {
         throw std::logic_error("Python object is not a dictionary. Cannot append");
      }

      PyObject *pKey, *pValue;

      for(auto mapIt = map.cbegin(); mapIt != map.cend(); ++mapIt)
      {
         pKey = PyUnicode_FromString(mapIt->second->tag().c_str());
         pValue = PyFloat_FromDouble(mapIt->second->value());
         PyDict_SetItem(pDict, pKey, pValue);
         Py_DECREF(pValue);
         Py_DECREF(pKey);
      }
   }

   void Tools::addDict(const std::map<std::string, MHDFloat>& map, PyObject* pDict)
   {
      if(!PyDict_Check(pDict))
      {
         throw std::logic_error("Python object is not a dictionary. Cannot append");
      }

      PyObject *pKey, *pValue;

      for(auto mapIt = map.cbegin(); mapIt != map.cend(); ++mapIt)
      {
         pKey = PyUnicode_FromString(mapIt->first.c_str());
         pValue = PyFloat_FromDouble(mapIt->second);
         PyDict_SetItem(pDict, pKey, pValue);
         Py_DECREF(pValue);
         Py_DECREF(pKey);
      }
   }

   void Tools::addDict(const std::map<std::string, std::string>& map, PyObject *pDict)
   {
      if(!PyDict_Check(pDict))
      {
         throw std::logic_error("Python object is not a dictionary. Cannot append");
      }

      PyObject *pKey, *pValue;

      for(auto mapIt = map.cbegin(); mapIt != map.cend(); ++mapIt)
      {
         pKey = PyUnicode_FromString(mapIt->first.c_str());
         pValue = PyUnicode_FromString(mapIt->second.c_str());
         PyDict_SetItem(pDict, pKey, pValue);
         Py_DECREF(pValue);
         Py_DECREF(pKey);
      }
   }

   void Tools::getList(std::vector<bool> &rList, PyObject *pList)
   {
      PyObject *pValue;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         bool isPeriodic = PyObject_IsTrue(pValue);

         rList.push_back(isPeriodic);
      }
   }

   void Tools::getList(std::vector<std::size_t> &rList, PyObject *pList)
   {
      PyObject *pValue, *pTmp;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         pTmp = PyUnicode_AsASCIIString(pValue);
         std::size_t id = Hasher::makeId(std::string(PyBytes_AsString(pTmp)));
         Py_DECREF(pTmp);

         rList.push_back(id);
      }
   }

   void Tools::getList(std::vector<std::string> &rList, PyObject *pList)
   {
      PyObject *pValue, *pTmp;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         pTmp = PyUnicode_AsASCIIString(pValue);
         std::string str = std::string(PyBytes_AsString(pTmp));
         Py_DECREF(pTmp);

         rList.push_back(str);
      }
   }

   void Tools::getList(std::vector<std::pair<std::size_t,FieldComponents::Spectral::Id> > &rList, PyObject *pList)
   {
      PyObject *pValue, *pTmp, *pTmp2;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         pTmp = PyTuple_GetItem(pValue,0);
         pTmp2 = PyUnicode_AsASCIIString(pTmp);
         std::size_t phys = Hasher::makeId(std::string(PyBytes_AsString(pTmp2)));
         Py_DECREF(pTmp2);

         pTmp = PyTuple_GetItem(pValue,1);
         pTmp2 = PyUnicode_AsASCIIString(pTmp);
         FieldComponents::Spectral::Id comp = QuICC::Tools::HumanToId::toComp(std::string(PyBytes_AsString(pTmp2)));
         Py_DECREF(pTmp2);

         rList.push_back(std::make_pair(phys, comp));
      }
   }

   void Tools::getDict(std::map<std::size_t,NonDimensional::SharedINumber> &rMap, PyObject *pDict, const bool replace)
   {
      PyObject *pValue, *pKey, *pList, *pTmp;

      pList = PyDict_Keys(pDict);

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pKey = PyList_GetItem(pList, i);
         pValue = PyDict_GetItem(pDict, pKey);

         pTmp = PyUnicode_AsASCIIString(pKey);
         std::size_t nd = Hasher::makeId(std::string(PyBytes_AsString(pTmp)));
         Py_DECREF(pTmp);

         if(replace && rMap.count(nd) > 0)
         {
            rMap.erase(nd);

         } else if(rMap.count(nd) > 0)
         {
            throw std::logic_error("Map key already existed!");
         }

         rMap.insert(std::make_pair(nd,NonDimensional::Coordinator::map().find(nd)->second->create(PyFloat_AsDouble(pValue))));
      }
   }

   void Tools::getDict(std::map<std::string,MHDFloat> &rMap, PyObject *pDict, const bool replace)
   {
      PyObject *pValue, *pKey, *pList, *pTmp;

      pList = PyDict_Keys(pDict);

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pKey = PyList_GetItem(pList, i);
         pValue = PyDict_GetItem(pDict, pKey);

         pTmp = PyUnicode_AsASCIIString(pKey);
         std::string str = std::string(PyBytes_AsString(pTmp));
         Py_DECREF(pTmp);

         if(replace && rMap.count(str) > 0)
         {
            rMap.find(str)->second = PyFloat_AsDouble(pValue);

         } else if(rMap.count(str) > 0)
         {
            throw std::logic_error("Map key already existed!");
         } else
         {
            rMap.insert(std::make_pair(str,PyFloat_AsDouble(pValue)));
         }
      }
   }

   PyObject* Tools::nparr2List(PyObject *pArr)
   {
      PyObject *pTmp, *pList;
      pTmp = PyObject_CallMethod(pArr, (char *)"flatten", "s", (char *)"F");
      pList = PyObject_CallMethod(pTmp, (char *)"tolist", NULL);
      Py_DECREF(pTmp);

      return pList;
   }

   PyObject* Tools::sparse2triplets(PyObject *pSpMat)
   {
      // Make sure matrix is in COO format
      PyObject *pTmp, *pArgs;
      PyObject *pCooMat;
      pTmp = PyObject_GetAttrString(pSpMat, (char *)"format");
      std::string format = std::string(PyBytes_AsString(PyUnicode_AsASCIIString(pTmp)));
      if(format == "coo")
      {
         pCooMat= pSpMat;
      }
      else
      {
         pCooMat = PyObject_CallMethod(pSpMat, (char *)"tocoo", NULL);
      }

      // Extra position of nonzero entries
      pArgs = PyTuple_New(3);
      pTmp = PyObject_GetAttrString(pCooMat, (char *)"row");
      PyTuple_SetItem(pArgs, 0, pTmp);
      pTmp = PyObject_GetAttrString(pCooMat, (char *)"col");
      PyTuple_SetItem(pArgs, 1, pTmp);
      pTmp = PyObject_GetAttrString(pCooMat, (char *)"data");
      PyTuple_SetItem(pArgs, 2, pTmp);
      PyObject *pBuiltins = PyEval_GetBuiltins(); 
      PyObject *pFct = PyDict_GetItemString(pBuiltins , "zip");
      if(! (pFct && PyCallable_Check(pFct)))
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
         throw std::logic_error("Python function loading error!");
      }
      PyObject *pZip = PyObject_CallObject(pFct, pArgs);
      Py_DECREF(pArgs);

      // Create list of triplets
      PyObject *pIterator = PyObject_GetIter(pZip);
      PyObject *pTriplet;
      PyObject *pList = PyList_New(0);
      while ((pTriplet = PyIter_Next(pIterator)))
      {
         PyList_Append(pList, pTriplet);
         Py_DECREF(pTriplet);
      }

      Py_DECREF(pIterator);

      return pList;
   }

   void Tools::fillArray(Array& rArray, PyObject* pPyArray)
   {
      PyObject *pValue;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyArray, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      Py_DECREF(pValue);

      // Allocate matrix
      rArray.resize(rows);

      // Convert Python matrix into a list
      pValue = Tools::nparr2List(pPyArray);

      long int count = 0;
      for (int i = 0; i < rows; i++)
      {
         rArray(i) = PyFloat_AsDouble(PyList_GetItem(pValue, count));
         count += 1;
      }

      Py_DECREF(pValue);
   }

   void Tools::fillMatrix(Matrix& rMatrix, PyObject* pPyMat)
   {
      PyObject *pValue;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));
      Py_DECREF(pValue);

      // Allocate matrix
      rMatrix.resize(rows,cols);

      // Convert Python matrix into a list
      pValue = Tools::nparr2List(pPyMat);

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

   void Tools::fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat)
   {
      PyObject *pValue, *pTmp;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));
      Py_DECREF(pValue);

      // Convert Python matrix into triplets
      pValue = Tools::sparse2triplets(pPyMat);

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

   void Tools::fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat)
   {
      PyObject *pValue, *pTmp;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));
      Py_DECREF(pValue);

      // Convert Python matrix into triplets
      pValue = Tools::sparse2triplets(pPyMat);

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

}
}
