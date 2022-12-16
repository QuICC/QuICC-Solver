/** 
 * @file PyModelBackend.cpp
 * @brief Source of the interface for Python model backend
 */

// First includes
//
#include "QuICC/PyQuICC/SystemHeader.hpp"

// Configuration includes
//

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Model/PyModelBackend.hpp"

// Project includes
//
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/PyQuICC/Tools.hpp"
#include "QuICC/PhysicalNames/Coordinator.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/Insulating.hpp"

namespace QuICC {

namespace Model {

   PyModelBackend::PyModelBackend(const std::string pyModule, const std::string pyClass)
   {
      // Initialize the python model wrapper
      this->mpWrapper = std::make_shared<PyQuICC::ModelWrapper>();
      this->mpWrapper->import(pyModule);
      this->mpWrapper->createModel(pyClass);
   }

   void PyModelBackend::enableGalerkin(const bool flag)
   {
      this->mpWrapper->enableGalerkin(flag);
   }

   std::vector<std::string> PyModelBackend::fieldNames() const
   {
      // Prepare Python call arguments
      PyObject *pValue;

      // Call model operator Python routine
      this->mpWrapper->setMethod((char *)"config_fields");
      pValue = this->mpWrapper->callMethod();

      // Create storage
      std::vector<std::string> names;
      PyQuICC::Tools::getList(names, pValue);

      // Clenup Python interpreter
      this->mpWrapper->cleanup();

      // Validate names
      this->checkFieldNames(names);

      return names;
   }

   std::vector<std::string> PyModelBackend::paramNames() const
   {
      // Prepare Python call arguments
      PyObject *pValue;

      // Call model operator Python routine
      this->mpWrapper->setMethod((char *)"nondimensional_parameters");
      pValue = this->mpWrapper->callMethod();

      // Create storage
      std::vector<std::string> names;
      PyQuICC::Tools::getList(names, pValue);

      // Cleanup Python interpreter
      this->mpWrapper->cleanup();

      // Validate names
      this->checkParamNames(names);

      return names;
   }

   std::vector<bool> PyModelBackend::isPeriodicBox() const
   {
      // Prepare Python call arguments
      PyObject *pValue;

      // Call model operator Python routine
      this->mpWrapper->setMethod((char *)"periodicity");
      pValue = this->mpWrapper->callMethod();

      // Create storage
      std::vector<bool> box;
      PyQuICC::Tools::getList(box, pValue);

      // Cleanup Python interpreter
      this->mpWrapper->cleanup();

      return box;
   }

   std::map<std::string,MHDFloat> PyModelBackend::automaticParameters(const std::map<std::string,MHDFloat>& cfg) const
   {
      //
      // Extend/modify parameters with automatically computed values
      //
      PyObject *pValue;
      PyObject *pTmp = PyQuICC::Tools::makeDict(cfg);

      // Call model operator Python routine
      PyObject *pArgs = PyTuple_New(1);
      PyTuple_SetItem(pArgs, 0, pTmp);
      this->mpWrapper->setMethod((char *)"automatic_parameters");
      pValue = this->mpWrapper->callMethod(pArgs);

      // Create storage
      std::map<std::string,MHDFloat> extra;
      PyQuICC::Tools::getDict(extra, pValue, true);
      Py_DECREF(pValue);
      Py_DECREF(pTmp);
      Py_DECREF(pArgs);

      // Cleanup Python interpreter
      this->mpWrapper->cleanup();

      return extra;
   }

   void PyModelBackend::equationInfo(bool& isComplex, SpectralFieldIds& im, SpectralFieldIds& exL, SpectralFieldIds& exNL, SpectralFieldIds& exNS, int& indexMode, const SpectralFieldId& fId, const Resolution& res) const
   {
      // Prepare Python call arguments
      PyObject *pArgs, *pTmp, *pValue;
      pArgs = PyTuple_New(2);

      // Get resolution
      pValue = PyQuICC::Tools::makeTuple(res.sim().dimensions(Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get field
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(PhysicalNames::Coordinator::tag(fId.first).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(QuICC::Tools::IdToHuman::toTag(fId.second).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 1, pTmp);

      // Call model operator Python routine
      this->mpWrapper->setMethod((char *)"equation_info");
      pValue = this->mpWrapper->callMethod(pArgs);
      Py_DECREF(pArgs);

      // Get Complex solver flag
      pTmp = PyTuple_GetItem(pValue, 0);
      isComplex = PyObject_IsTrue(pTmp);

      // Get Implicit fields
      pTmp = PyTuple_GetItem(pValue, 1);
      PyQuICC::Tools::getList(im, pTmp);

      // Get Explicit linear fields
      pTmp = PyTuple_GetItem(pValue, 2);
      PyQuICC::Tools::getList(exL, pTmp);

      // Get Explicit nonlinear fields
      pTmp = PyTuple_GetItem(pValue, 3);
      PyQuICC::Tools::getList(exNL, pTmp);

      // Get Explicit nextstep fields
      pTmp = PyTuple_GetItem(pValue, 4);
      PyQuICC::Tools::getList(exNS, pTmp);

      // Get index mode
      pTmp = PyTuple_GetItem(pValue, 5);
      indexMode = PyLong_AsLong(pTmp);

      // Finalise Python interpreter
      Py_DECREF(pValue);
   }

   std::map<std::string,int> PyModelBackend::getPyBcMap(const BcMap& bcs) const
   {
      std::map<std::string,int> m;

      for(auto bc: bcs)
      {
         auto physId = PhysicalNames::Coordinator::tag(bc.first);
         int bcVal = -1;
         if(bc.second == Bc::Name::NoSlip::id())
         {
            bcVal = 0;
         }
         else if(bc.second == Bc::Name::StressFree::id())
         {
            bcVal = 1;
         }
         else if(bc.second == Bc::Name::FixedTemperature::id())
         {
            bcVal = 0;
         }
         else if(bc.second == Bc::Name::FixedFlux::id())
         {
            bcVal = 1;
         }
         else if(bc.second == Bc::Name::Insulating::id())
         {
            bcVal = 0;
         }
         else
         {
            throw std::logic_error("Conversion of boundary condition name to python flag is unknown");
         }

         m.emplace(std::pair(physId, bcVal));
      }

      return m;
   }

   void PyModelBackend::operatorInfo(ArrayI& tauN, ArrayI& galN, MatrixI& galShift, ArrayI& rhsCols, ArrayI& sysN, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const
   {
      PyObject *pArgs, *pTmp, *pValue;

      // Prepare Python call
      pArgs = PyTuple_New(4);

      // Get boundary conditions
      auto pyBcs = this->getPyBcMap(bcs);
      pValue = PyQuICC::Tools::makeDict(pyBcs);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get field
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(PhysicalNames::Coordinator::tag(fId.first).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(QuICC::Tools::IdToHuman::toTag(fId.second).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 3, pTmp);

      // Call model operator Python routine
      this->mpWrapper->setMethod((char *)"operator_info");

      // Loop overall matrices/eigs
      for(int matIdx = 0; matIdx < tauN.size(); ++matIdx)
      {
         // Get resolution
         auto eigs = coupling.getIndexes(res, matIdx);
         pValue = PyQuICC::Tools::makeTuple(res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0)));
         PyTuple_SetItem(pArgs, 0, pValue);

         // Get the eigen direction values
         pValue = PyQuICC::Tools::makeTuple(eigs);
         PyTuple_SetItem(pArgs, 1, pValue);

         pValue = this->mpWrapper->callMethod(pArgs);

         // Get block information
         pTmp = PyTuple_GetItem(pValue, 0);
         tauN(matIdx) = PyLong_AsLong(pTmp);
         pTmp = PyTuple_GetItem(pValue, 1);
         galN(matIdx) = PyLong_AsLong(pTmp);
         pTmp = PyTuple_GetItem(pValue, 2);
         galShift(matIdx,0) = PyLong_AsLong(PyTuple_GetItem(pTmp, 0));
         galShift(matIdx,1) = PyLong_AsLong(PyTuple_GetItem(pTmp, 1));
         galShift(matIdx,2) = PyLong_AsLong(PyTuple_GetItem(pTmp, 2));
         pTmp = PyTuple_GetItem(pValue, 3);
         rhsCols(matIdx) = PyLong_AsLong(pTmp);
         pTmp = PyTuple_GetItem(pValue, 4);
         sysN(matIdx) = PyLong_AsLong(pTmp);

         Py_DECREF(pValue);
      }
      Py_DECREF(pArgs);
   }

   PyObject*  PyModelBackend::baseArguments(const int tupleSize, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, BcMap bcs, const NonDimensional::NdMap& nds) const
   {
      // Prepare Python call arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(tupleSize);

      // Get resolution
      pValue = PyQuICC::Tools::makeTuple(res.counter().dimensions(Dimensions::Space::SPECTRAL, eigs.at(0)));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Get equation parameters
      pValue = PyQuICC::Tools::makeDict(nds);
      PyTuple_SetItem(pArgs, 1, pValue);

      // Get the eigen direction values
      pValue = PyQuICC::Tools::makeTuple(eigs);
      PyTuple_SetItem(pArgs, 2, pValue);

      // Get boundary conditions
      auto pyBcs = this->getPyBcMap(bcs);
      pValue = PyQuICC::Tools::makeDict(pyBcs);
      // ... append boundary type
      std::map<std::string,std::string> bcsStr;
      bcsStr.insert(std::make_pair("bcType", ModelOperatorBoundary::Coordinator::tag(bcType)));
      PyQuICC::Tools::addDict(bcsStr, pValue);

      PyTuple_SetItem(pArgs, 3, pValue);

      return pArgs;
   }

   void PyModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      // Get first four standard arguments in a tuple of size 6
      PyObject *pArgs;
      int argStart;
      pArgs = this->baseArguments(6, bcType, res, eigs, bcs, nds);
      argStart = 4;

      // Prepare Python call arguments
      PyObject *pTmp, *pValue, *pList;

      // Get list of implicit fields
      pList = PyList_New(0);
      for(auto fIt = imRange.first; fIt != imRange.second; ++fIt)
      {
         pTmp = PyTuple_New(2);
         pValue = PyUnicode_FromString(PhysicalNames::Coordinator::tag(fIt->first).c_str());
         PyTuple_SetItem(pTmp, 0, pValue);
         pValue = PyUnicode_FromString(QuICC::Tools::IdToHuman::toTag(fIt->second).c_str());
         PyTuple_SetItem(pTmp, 1, pValue);

         PyList_Append(pList, pTmp);
         Py_DECREF(pTmp);
      }
      PyTuple_SetItem(pArgs, argStart, pList);

      // Set the restriction option
      #ifdef QUICC_MPI
         #ifdef QUICC_MPISPSOLVE
            std::vector<int> slow;
            std::vector<std::vector<int> > middle;

            res.buildRestriction(slow, middle, matIdx);
            PyObject *pSlow, *pMiddle;

            if(middle.size() > 0)
            {
               pSlow = PyQuICC::Tools::makeList(slow);
               pMiddle = PyQuICC::Tools::makeList(middle);
               pTmp = PyTuple_New(2);
               PyTuple_SetItem(pTmp, 0, pSlow);
               PyTuple_SetItem(pTmp, 1, pMiddle);
            } else
            {
               pTmp = PyQuICC::Tools::makeList(slow);
            }

            PyTuple_SetItem(pArgs, argStart+1, pTmp);
         #else
            // MPI code with serial sparse solver
            Py_INCREF(Py_None);
            PyTuple_SetItem(pArgs, argStart+1, Py_None);
         #endif //QUICC_MPISPSOLVE
      #else
         // Serial code can't have any restriction
         Py_INCREF(Py_None);
         PyTuple_SetItem(pArgs, argStart+1, Py_None);
      #endif //QUICC_MPI

      // Call model operator Python routine
      this->mpWrapper->setMethod(QuICC::ModelOperator::Coordinator::tag(opId));
      pValue = this->mpWrapper->callMethod(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      this->mpWrapper->fillMatrix(rModelMatrix, pValue);
      Py_DECREF(pValue);

      // Finalise Python interpreter
      this->mpWrapper->ModelWrapper::cleanup();
   }

   void PyModelBackend::galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      // Get first four standard arguments in a tuple of size 5
      PyObject *pArgs = this->baseArguments(6, ModelOperatorBoundary::Stencil::id(), res, eigs, bcs, nds);

      // Prepare Python call arguments
      PyObject *pTmp, *pValue;

      // Get field
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(PhysicalNames::Coordinator::tag(fId.first).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(QuICC::Tools::IdToHuman::toTag(fId.second).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 4, pTmp);
      if(makeSquare)
      {
         Py_INCREF(Py_True);
         PyTuple_SetItem(pArgs, 5, Py_True);
      } else
      {
         Py_INCREF(Py_False);
         PyTuple_SetItem(pArgs, 5, Py_False);
      }

      // Call model operator Python routine
      this->mpWrapper->setMethod(ModelOperator::Coordinator::tag(ModelOperator::Stencil::id()));
      pValue = this->mpWrapper->callMethod(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      this->mpWrapper->fillMatrix(mat, pValue);
      Py_DECREF(pValue);

      // Finalise Python interpreter
      this->mpWrapper->cleanup();
   }

   void PyModelBackend::explicitBlock(DecoupledZSparse& mat, const SpectralFieldId& fId, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      // Get first four standard arguments in a tuple of size 7
      PyObject *pArgs = this->baseArguments(7, ModelOperatorBoundary::FieldToRhs::id(), res, eigs, bcs, nds);

      // Prepare Python call arguments
      PyObject *pTmp, *pValue;

      // Get field row
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(PhysicalNames::Coordinator::tag(fId.first).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(QuICC::Tools::IdToHuman::toTag(fId.second).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 4, pTmp);

      // Get field col
      pTmp = PyTuple_New(2);
      pValue = PyUnicode_FromString(PhysicalNames::Coordinator::tag(fieldId.first).c_str());
      PyTuple_SetItem(pTmp, 0, pValue);
      pValue = PyUnicode_FromString(QuICC::Tools::IdToHuman::toTag(fieldId.second).c_str());
      PyTuple_SetItem(pTmp, 1, pValue);
      PyTuple_SetItem(pArgs, 5, pTmp);

      // Set the restriction option
      #ifdef QUICC_MPI
         #ifdef QUICC_MPISPSOLVE
            std::vector<int> slow;
            std::vector<std::vector<int> > middle;

            res.buildRestriction(slow, middle, matIdx);
            PyObject *pSlow, *pMiddle;

            if(middle.size() > 0)
            {
               pSlow = PyQuICC::Tools::makeList(slow);
               pMiddle = PyQuICC::Tools::makeList(middle);
               pTmp = PyTuple_New(2);
               PyTuple_SetItem(pTmp, 0, pSlow);
               PyTuple_SetItem(pTmp, 1, pMiddle);
            } else
            {
               pTmp = PyQuICC::Tools::makeList(slow);
            }

            PyTuple_SetItem(pArgs, 6, pTmp);
         #else
            // MPI code with serial sparse solver
            Py_INCREF(Py_None);
            PyTuple_SetItem(pArgs, 6, Py_None);
         #endif //QUICC_MPISPSOLVE
      #else
         // Serial code can't have any restriction
         Py_INCREF(Py_None);
         PyTuple_SetItem(pArgs, 6, Py_None);
      #endif //QUICC_MPI

      // Call model operator Python routine
      this->mpWrapper->setMethod(QuICC::ModelOperator::Coordinator::tag(opId));
      pValue = this->mpWrapper->callMethod(pArgs);
      Py_DECREF(pArgs);

      // Convert Python matrix into triplets
      this->mpWrapper->fillMatrix(mat, pValue);
      Py_DECREF(pValue);

      // Finalise Python interpreter
      this->mpWrapper->cleanup();
   }

}
}
