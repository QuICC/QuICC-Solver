/**
 * @file Cartesian1DNusseltDZWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer through the Z boundary extracted from temperature field
 */

// Configuration includes
//

// System includes
//
#include <iomanip>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/Cartesian1DNusseltDZWriter.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/NusseltTags.hpp"
#include "QuICC/PyQuICC/CoreWrapper.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   Cartesian1DNusseltDZWriter::Cartesian1DNusseltDZWriter(std::string type)
      : IVariableAsciiWriter(Tags::Nusselt::BASENAME, Tags::Nusselt::EXTENSION, Tags::Nusselt::HEADER, type, Tags::Nusselt::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   Cartesian1DNusseltDZWriter::~Cartesian1DNusseltDZWriter()
   {
   }

   void Cartesian1DNusseltDZWriter::init()
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      scalar_iterator sit = sRange.first;

      // Initialise python wrapper
      PyQuICC::CoreWrapper::init();
      PyQuICC::CoreWrapper::import("quicc.geometry.cartesian.cartesian_1d");

      // Prepare arguments
      PyObject *pArgs, *pValue, *pTmp;
      pArgs = PyTuple_New(4);

      // Get resolution
      pValue = PyLong_FromLong(sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);
      pTmp = PyFloat_FromDouble(this->mPhysical.find(NonDimensional::Lower1d::id())->second->value());
      PyTuple_SetItem(pArgs, 1, pTmp);
      pTmp = PyFloat_FromDouble(this->mPhysical.find(NonDimensional::Upper1d::id())->second->value());
      PyTuple_SetItem(pArgs, 2, pTmp);

      // Set boundary condition
      pValue = PyDict_New();
      PyDict_SetItem(pValue, PyLong_FromLong(0), PyLong_FromLong(21));
      Py_INCREF(Py_False);
      PyDict_SetItem(pValue, PyUnicode_FromString("use_parity"), Py_False);
      PyTuple_SetItem(pArgs, 3, pValue);

      // Call zblk and use derivative boundary condition
      PyQuICC::CoreWrapper::setFunction("zblk");
      pValue = PyQuICC::CoreWrapper::callFunction(pArgs);

      // Fill matrix and cleanup
      PyQuICC::CoreWrapper::fillMatrix(this->mNusseltOp, pValue);
      Py_DECREF(pValue);
      PyQuICC::CoreWrapper::finalize();

      IVariableAsciiWriter::init();
   }

   void Cartesian1DNusseltDZWriter::write()
   {
      // Create file
      this->preWrite();

      Cartesian1DNusseltDZWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      Cartesian1DNusseltDZWriter::scalar_iterator sit = sRange.first;

      ArrayI mode = sit->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->mode(0);
      MHDFloat nusselt = 0.0;
      if(mode(2) == 0 && mode(3) == 0)
      {
         // Compute Nusselt number
         nusselt = -(this->mNusseltOp*sit->second->dom(0).perturbation().profile(0,0).real()).row(0)(0);
      }

      // Get the "global" Nusselt number from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &nusselt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << 1.0 + nusselt << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if Nusselt number is NaN
      if(std::isnan(nusselt))
      {
         QuICCEnv().abort("Nusselt number is NaN!");
      }
   }

}
}
}
