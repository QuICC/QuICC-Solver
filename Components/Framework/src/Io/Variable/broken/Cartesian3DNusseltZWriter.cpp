/**
 * @file Cartesian3DNusseltZWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer for a 3D box throught Z
 */

// Configuration includes
//

// System includes
//
#include <iomanip>

// External includes
//

// Class include
//
#include "QuICC/Io/Variable/Cartesian3DNusseltZWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/NusseltTags.hpp"
#include "QuICC/PyQuICC/CoreWrapper.hpp"
#include "QuICC/DecoupledComplexInternal.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   Cartesian3DNusseltZWriter::Cartesian3DNusseltZWriter(std::string type)
      : IVariableAsciiWriter(Tags::Nusselt::BASENAME, Tags::Nusselt::EXTENSION, Tags::Nusselt::HEADER, type, Tags::Nusselt::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   Cartesian3DNusseltZWriter::~Cartesian3DNusseltZWriter()
   {
   }

   void Cartesian3DNusseltZWriter::init()
   {
      Cartesian3DNusseltZWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      Cartesian3DNusseltZWriter::scalar_iterator sit = sRange.first;

      // Initialise python wrapper
      PyQuICC::CoreWrapper::init();
      PyQuICC::CoreWrapper::import("quicc.geometry.cartesian.cartesian_3d");

      // Prepare arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(4);

      // Get resolution
      pValue = PyLong_FromLong(sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);
      pValue = PyLong_FromLong(sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 1, pValue);
      pValue = PyLong_FromLong(sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 2, pValue);
      // Set scale
      pValue = PyLong_FromLong(2.0);
      PyTuple_SetItem(pArgs, 3, pValue);

      // Call avgFlux_z
      PyQuICC::CoreWrapper::setFunction("avgFlux_z");
      pValue = PyQuICC::CoreWrapper::callFunction(pArgs);

      // Fill matrix and clenup
      PyQuICC::CoreWrapper::fillMatrix(this->mNusseltOp, pValue);
      Py_DECREF(pValue);
      PyQuICC::CoreWrapper::finalize();

      IVariableAsciiWriter::init();
   }

   void Cartesian3DNusseltZWriter::write()
   {
      // Create file
      this->preWrite();

      Cartesian3DNusseltZWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      Cartesian3DNusseltZWriter::scalar_iterator sit = sRange.first;

      // Copy data
      int l;
      int k_;
      int j_;
      int dimK = sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL)*sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
      int dimJ = sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      Array field(dimK*sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL));
      for(int k = 0; k < sit->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); k++)
      {
         k_ = sit->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)*dimK;
         for(int j = 0; j < sit->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
         {
            j_ = sit->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k)*dimJ;
            for(int i = 0; i < sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL); i++)
            {
               // Compute correct position
               l = k_ + j_ + i;

               // Copy field value into field
               Datatypes::internal::setScalar(field, l, sit->second->dom(0).perturbation().point(i,j,k));
            }
         }
      }

      // Get the "global" field
      #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
         MPI_Reduce(MPI_IN_PLACE, field.data(), field.rows(), MPI_DOUBLE, MPI_SUM, 0, QuICC::FramworkMacro::getSubComm(MpiFramwork::LOCAL));
      #endif //defined QUICC_MPI && defined QUICC_MPISPSOLVE

      Matrix nusselt = -this->mNusseltOp*field;
      assert(nusselt.rows() == nusselt.cols() && nusselt.rows() == 1);
      nusselt(0,0) += 1;

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << nusselt(0,0) << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
}
