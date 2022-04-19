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

}
}
