/**
 * @file ToolsTest.cpp
 * @brief Tests for the Tools wrapper
 */


// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/PyQuICC/Coordinator.hpp"
#include "QuICC/PyQuICC/Tools.hpp"
#include "QuICC/PyQuICC/CoreWrapper.hpp"

namespace currentts = QuICC::PyQuICC;

TEST_CASE( "PyQuICC tools", "[Tools]" ){

   // Initialize Python interpreter
   currentts::Coordinator::init();

   std::string mod = "quicc.testsuite.pyquicc.tools_tests";
   currentts::CoreWrapper::import(mod);

   //
   {
      INFO( "Check convertion to tuple" );
      QuICC::ArrayI a(4);
      a << 1, 2, 3, 5;
      PyObject *p1, *p2;
      p1 = currentts::Tools::makeTuple(a);
      std::string fct = "check_tuple";
      currentts::CoreWrapper::setFunction(fct);
      PyObject *pArgs = PyTuple_New(1);
      Py_INCREF(p1);
      PyTuple_SetItem(pArgs, 0, p1);
      p2 = currentts::CoreWrapper::callFunction(pArgs);
      CHECK(PyObject_IsTrue(p2));
      Py_DECREF(pArgs);
      Py_DECREF(p1);
      Py_DECREF(p2);
   }

   //
   {
      INFO( "Check NumPy array to list" );
      std::string fct = "arr";
      PyObject *p1, *p2, *pArgs;
      currentts::CoreWrapper::setFunction(fct);
      p1 = currentts::CoreWrapper::callFunction();
      p2 = currentts::Tools::nparr2List(p1);
      Py_DECREF(p1);
      fct = "check_arr2list";
      currentts::CoreWrapper::setFunction(fct);
      pArgs = PyTuple_New(1);
      Py_INCREF(p2);
      PyTuple_SetItem(pArgs, 0, p2);
      p1 = currentts::CoreWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);
      Py_DECREF(p2);
      CHECK(PyObject_IsTrue(p1));
      Py_DECREF(p1);
   }

   //
   {
      INFO( "Check NumPy array to list" );
      std::string fct = "mat";
      PyObject *p1, *p2, *pArgs;
      currentts::CoreWrapper::setFunction(fct);
      p1 = currentts::CoreWrapper::callFunction();
      p2 = currentts::Tools::nparr2List(p1);
      Py_DECREF(p1);
      fct = "check_mat2list";
      currentts::CoreWrapper::setFunction(fct);
      pArgs = PyTuple_New(1);
      Py_INCREF(p2);
      PyTuple_SetItem(pArgs, 0, p2);
      p1 = currentts::CoreWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);
      Py_DECREF(p2);
      CHECK( PyObject_IsTrue(p1) );
      Py_DECREF(p1);
   }

   //
   {
      INFO( "Convert NumPy array to QuICC Array" );
      std::string fct = "arr";
      PyObject *p1;
      currentts::CoreWrapper::setFunction(fct);
      p1 = currentts::CoreWrapper::callFunction();
      QuICC::Array qArr;
      currentts::Tools::fillArray(qArr, p1);
      Py_DECREF(p1);
      QuICC::Array qArrRef(4);
      qArrRef << 0.1, 0.2, 0.3, 0.5;
      bool sameSize = (qArr.size() == qArrRef.size());
      CHECK( sameSize );
      if(sameSize)
      {
         CHECK( (qArr == qArrRef) );
      }
   }

   //
   {
      INFO( "Convert NumPy 2D array to QuICC Matrix" )
      std::string fct = "mat";
      PyObject *p1;
      currentts::CoreWrapper::setFunction(fct);
      p1 = currentts::CoreWrapper::callFunction();
      QuICC::Matrix qMat;
      currentts::Tools::fillMatrix(qMat, p1);
      Py_DECREF(p1);
      QuICC::Matrix qMatRef(3,4);
      qMatRef << 0.1, 0.2, 0.3, 0.5,
              1.1, 1.2, 1.3, 1.5,
              2.1, 2.2, 2.3, 2.5;
      bool sameSize = (qMat.size() == qMatRef.size());
      CHECK( sameSize );
      CHECK( (qMat == qMatRef) );
   }

   //
   {
      INFO( "Check SciPy sparse matrix to triplets" );
      std::string fct = "spMat";
      PyObject *p1, *p2, *pArgs;
      currentts::CoreWrapper::setFunction(fct);
      p1 = currentts::CoreWrapper::callFunction();
      p2 = currentts::Tools::sparse2triplets(p1);
      Py_DECREF(p1);
      fct = "check_sparse2triplets";
      currentts::CoreWrapper::setFunction(fct);
      pArgs = PyTuple_New(1);
      Py_INCREF(p2);
      PyTuple_SetItem(pArgs, 0, p2);
      p1 = currentts::CoreWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);
      Py_DECREF(p2);
      CHECK(PyObject_IsTrue(p1));
      Py_DECREF(p1);
   }

   //
   {
      INFO( "Check SciPy sparse matrix (not COO) to triplets" );
      std::string fct = "spMat2";
      PyObject *p1, *p2, *pArgs;
      currentts::CoreWrapper::setFunction(fct);
      p1 = currentts::CoreWrapper::callFunction();
      p2 = currentts::Tools::sparse2triplets(p1);
      Py_DECREF(p1);
      fct = "check_sparse2triplets";
      currentts::CoreWrapper::setFunction(fct);
      pArgs = PyTuple_New(1);
      Py_INCREF(p2);
      PyTuple_SetItem(pArgs, 0, p2);
      p1 = currentts::CoreWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);
      Py_DECREF(p2);
      CHECK(PyObject_IsTrue(p1));
      Py_DECREF(p1);
   }

   //
   {
      INFO( "Check SciPy sparse matrix to QuICC SparseMatrix" );
      std::string fct = "spMat";
      PyObject *p1;
      currentts::CoreWrapper::setFunction(fct);
      p1 = currentts::CoreWrapper::callFunction();
      QuICC::SparseMatrix qSpMat;
      currentts::Tools::fillMatrix(qSpMat, p1);
      Py_DECREF(p1);
      QuICC::SparseMatrix qSpMatRef(5,5);
      for(int i = 0; i < qSpMatRef.rows(); i++)
      {
         qSpMatRef.coeffRef(i,i) = static_cast<double>(i+1);
         if(i < qSpMatRef.rows()-1)
         {
            qSpMatRef.coeffRef(i,i+1) = -static_cast<double>(i+1);
         }
      }
      bool sameSize = (qSpMat.size() == qSpMatRef.size());
      CHECK( sameSize );
      for(int j = 0; j < qSpMatRef.cols(); j++)
      {
         for(int i = 0; i < qSpMatRef.rows(); i++)
         {
            CHECK( qSpMat.coeff(i,j) == qSpMatRef.coeff(i,j) );
         }
      }
   }

   //
   {
      INFO("Check SciPy sparse complex matrix to QuICC DecoupledZSparse" );
      std::string fct = "spMatZ";
      PyObject *p1;
      currentts::CoreWrapper::setFunction(fct);
      p1 = currentts::CoreWrapper::callFunction();
      QuICC::DecoupledZSparse qSpMatZ;
      currentts::Tools::fillMatrix(qSpMatZ, p1);
      Py_DECREF(p1);
      QuICC::DecoupledZSparse qSpMatRefZ(5,5);
      for(int i = 0; i < qSpMatRefZ.real().rows(); i++)
      {
         qSpMatRefZ.real().coeffRef(i,i) = static_cast<double>(i+1);
         qSpMatRefZ.imag().coeffRef(i,i) = static_cast<double>(i+11);
         if(i < qSpMatRefZ.real().rows()-1)
         {
            qSpMatRefZ.real().coeffRef(i,i+1) = -static_cast<double>(i+1);
            qSpMatRefZ.imag().coeffRef(i,i+1) = -static_cast<double>(i+7);
         }
      }
      bool sameSize = (qSpMatZ.real().size() == qSpMatRefZ.real().size());
      CHECK( sameSize );
      sameSize = (qSpMatZ.imag().size() == qSpMatRefZ.imag().size());
      CHECK( sameSize );
      for(int j = 0; j < qSpMatRefZ.real().cols(); j++)
      {
         for(int i = 0; i < qSpMatRefZ.real().rows(); i++)
         {
            CHECK( qSpMatZ.real().coeff(i,j) == qSpMatRefZ.real().coeff(i,j) );
            CHECK( qSpMatZ.imag().coeff(i,j) == qSpMatRefZ.imag().coeff(i,j) );
         }
      }
   }

   // Finalize Python interpreter
   currentts::Coordinator::finalize();
}
