/**
 * @file Cartesian1DNusseltXWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer through the X boundary
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
#include "QuICC/Io/Variable/Cartesian1DNusseltXWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/ScalarFields/FieldTools.hpp"
#include "QuICC/Io/Variable/NusseltTags.hpp"
#include "QuICC/PyQuICC/CoreWrapper.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   Cartesian1DNusseltXWriter::Cartesian1DNusseltXWriter(std::string type)
      : IVariableAsciiWriter(Tags::Nusselt::BASENAME, Tags::Nusselt::EXTENSION, Tags::Nusselt::HEADER, type, Tags::Nusselt::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   Cartesian1DNusseltXWriter::~Cartesian1DNusseltXWriter()
   {
   }

   void Cartesian1DNusseltXWriter::init()
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 2);
      scalar_iterator sit = sRange.first;

      // Initialise python wrapper
      PyQuICC::CoreWrapper::init();
      PyQuICC::CoreWrapper::import("quicc.geometry.cartesian.cartesian_1d");

      // Prepare arguments
      PyObject *pArgs, *pValue;
      pArgs = PyTuple_New(1);

      // Get resolution
      pValue = PyLong_FromLong(sit->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call avg
      PyQuICC::CoreWrapper::setFunction("avg");
      pValue = PyQuICC::CoreWrapper::callFunction(pArgs);

      // Fill matrix and cleanup
      PyQuICC::CoreWrapper::fillMatrix(this->mNusseltOp, pValue);
      Py_DECREF(pValue);
      PyQuICC::CoreWrapper::finalize();

      IVariableAsciiWriter::init();
   }

   void Cartesian1DNusseltXWriter::write()
   {
      // Create file
      this->preWrite();

      Cartesian1DNusseltXWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 2);
      Cartesian1DNusseltXWriter::scalar_iterator nlIt;
      Cartesian1DNusseltXWriter::scalar_iterator tIt;
      Cartesian1DNusseltXWriter::scalar_iterator it = sRange.first;
      if(it->first == PhysicalNames::Temperature::id())
      {
         tIt = it;
         ++it;
         nlIt = it;
      } else
      {
         nlIt = it;
         ++it;
         tIt = it;
      }

      ArrayI mode = tIt->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->mode(0);
      ArrayZ tmp;
      Array nusselt = Array::Zero(2);
      if(mode(2) == 0 && mode(3) == 0)
      {
         // Heat flux
         nusselt(0) = nlIt->second->dom(0).perturbation().point(0,0,0).real();

         // Mean heat
         tmp = this->mNusseltOp*tIt->second->dom(0).perturbation().profile(0,0);
         nusselt(1) = tmp(0).real();
      }

//      Array nusselt = Array::Zero(3);
//      ArrayZ tmp;
//      MHDComplex zk_;
//      for(int k = 0; k < tIt->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
//      {
//         if(tIt->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
//         {
//            nusselt(0) += nlIt->second->dom(0).perturbation().point(0,0,k).real();
//
//            if(tIt->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k) <= tIt->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)/2)
//            {
//               zk_ = MHDComplex(0.0, static_cast<MHDFloat>(tIt->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)));
//            } else
//            {
//               zk_ = MHDComplex(0.0, static_cast<MHDFloat>(tIt->second->dom(0).res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k)) - tIt->second->dom(0).res().sim().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL));
//            }
//            tmp = this->mNusseltOp*tIt->second->dom(0).perturbation().profile(0,k);
//            nusselt(1) += tmp(0).real();
//            nusselt(2) += (zk_*tmp(0)).real();
//         }
//      }

      // Get the "global" Nusselt number from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, nusselt.data(), nusselt.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << 1.0 - nusselt(0) << "\t" << nusselt(1) << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if Nusselt number is NaN
      if(std::isnan(nusselt(0)) || std::isnan(nusselt(1)))
      {
         QuICCEnv().abort(99);

         throw std::logic_error("Nusselt number is NaN!");
      }
   }

}
}
}
