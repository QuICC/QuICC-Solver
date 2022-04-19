/**
 * @file Coordinator.cpp
 * @brief Source of Python interpreter wrapper
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PyQuICC/Coordinator.hpp"

// Project includes
//
#include "QuICC/PyQuICC/Config.hpp"

namespace QuICC {

namespace PyQuICC {

   int Coordinator::sCounter = 0;

   void Coordinator::init()
   {
      if(Coordinator::sCounter == 0)
      {
         // Initialize the python interpreter
         Py_Initialize();

         // Setup the search path
         PyObject* sysPath = PySys_GetObject((char*)"path");
         PyList_Append(sysPath, PyUnicode_FromString(QUICC_PYTHON_DIR));
      }

      Coordinator::registerWrapper();
   }

   void Coordinator::registerWrapper()
   {
      Coordinator::sCounter++;
   }

   void Coordinator::unregisterWrapper()
   {
      Coordinator::sCounter--;
   }

   void Coordinator::finalize()
   {
      Coordinator::unregisterWrapper();

      if(Coordinator::sCounter == 0)
      {
         // Finalize
         Py_Finalize();
      }
   }

   Coordinator::Coordinator()
   {
   }

   Coordinator::~Coordinator()
   {
   }

}
}
