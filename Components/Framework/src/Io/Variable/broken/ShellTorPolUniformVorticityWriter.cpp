/**
 * @file ShellTorPolUniformVorticityWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics uniform vorticity calculation for scalar field in a spherical shell
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

// Configuration includes
//

// System includes
//
#include <iomanip>
#include <assert.h>
// External includes
//

// Class include
//
#include "QuICC/Io/Variable/ShellTorPolUniformVorticityWriter.hpp"

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Ro.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Io/Variable/AverageTags.hpp"
#include "QuICC/TypeSelectors/ScalarSelector.hpp"
#include "QuICC/PyQuICC/CoreWrapper.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

   ShellTorPolUniformVorticityWriter::ShellTorPolUniformVorticityWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + AverageTags::BASENAME, AverageTags::EXTENSION, prefix + AverageTags::HEADER, type, AverageTags::VERSION, Dimensions::Space::SPECTRAL, IVariableAsciiWriter::EXTEND)
   {
   }

   ShellTorPolUniformVorticityWriter::~ShellTorPolUniformVorticityWriter()
   {
   }

   void ShellTorPolUniformVorticityWriter::init()
   {
      // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
      MHDFloat ro = this->mPhysical.find(NonDimensional::Ro::id())->second->value();
      MHDFloat rratio = this->mPhysical.find(NonDimensional::RRatio::id())->second->value();
      MHDFloat ri = ro*rratio;

      this->mDelta = this->mPhysical.find(NonDimensional::Ekman::id())->second->value();
      this->mDelta = std::pow(this->mDelta,0.5)*10.;
      //this->mVolume = (std::pow(ro-this->mDelta,5)-std::pow(ri+this->mDelta,5))/(5.* std::sqrt(3.));
      this->mVolume = (std::pow(ro-this->mDelta,5)-std::pow(ri+this->mDelta,5))/(5.* std::sqrt(3./(4*Math::PI)));

      // Initialise python wrapper
      PyQuICC::CoreWrapper::init();
      PyQuICC::CoreWrapper::import("quicc.geometry.spherical.shell_radius");

      // Prepare arguments
      PyObject *pArgs, *pValue, *pBC;
      pArgs = PyTuple_New(4);

      // ... compute a, b factors
      PyObject *pTmp = PyTuple_New(2);
      PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(ro));
      PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(rratio));
      PyQuICC::CoreWrapper::setFunction("linear_r2x");
      pValue = PyQuICC::CoreWrapper::callFunction(pTmp);
      PyTuple_SetItem(pArgs, 1, PyTuple_GetItem(pValue, 0));
      PyTuple_SetItem(pArgs, 2, PyTuple_GetItem(pValue, 1));
      // ... create boundray condition (none)
      pBC = PyDict_New();
      PyDict_SetItem(pBC, PyLong_FromLong(0), PyLong_FromLong(0));
      PyDict_SetItem(pBC, PyUnicode_FromString("cr"), PyLong_FromLong(2));
      PyTuple_SetItem(pArgs, 3, pBC);

      // Get resolution
      int nR = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      pValue = PyLong_FromLong(nR+2);
      PyTuple_SetItem(pArgs, 0, pValue);

      // Call r^2
      PyQuICC::CoreWrapper::setFunction("r2");
      pValue = PyQuICC::CoreWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpR2(nR+2, nR);
      PyQuICC::CoreWrapper::fillMatrix(tmpR2, pValue);
      Py_DECREF(pValue);

      // change a couple of arguments
      pValue = PyLong_FromLong(nR+3);
      PyTuple_SetItem(pArgs, 0, pValue);
      pBC = PyDict_New();
      PyDict_SetItem(pBC, PyLong_FromLong(0), PyLong_FromLong(0));
      PyDict_SetItem(pBC, PyUnicode_FromString("cr"), PyLong_FromLong(1));
      PyTuple_SetItem(pArgs, 3, pBC);
      // Call r^1
      PyQuICC::CoreWrapper::setFunction("r1");
      pValue = PyQuICC::CoreWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpR1(nR+3, nR+2);
      PyQuICC::CoreWrapper::fillMatrix(tmpR1, pValue);
      Py_DECREF(pValue);

      // change a couple of arguments
      pValue = PyLong_FromLong(nR+4);
      PyTuple_SetItem(pArgs, 0, pValue);
      // Call i1
      PyQuICC::CoreWrapper::setFunction("i1");
      pValue = PyQuICC::CoreWrapper::callFunction(pArgs);
      // Fill matrix and cleanup
      SparseMatrix tmpI1(nR+4, nR+3);
      PyQuICC::CoreWrapper::fillMatrix(tmpI1, pValue);
      Py_DECREF(pValue);

      // prepare the points of cut offs
      PyObject * pList;
      pList = PyList_New(2);
      PyList_SetItem(pList, 0, PyFloat_FromDouble(ro-this->mDelta));
      PyList_SetItem(pList, 1, PyFloat_FromDouble(ri+this->mDelta));
      PyTuple_SetItem(pArgs, 3, pList);

      //Call proj
      PyQuICC::CoreWrapper::import("quicc.projection.shell");
      PyQuICC::CoreWrapper::setFunction("proj");
      pValue = PyQuICC::CoreWrapper::callFunction(pArgs);
      SparseMatrix tmpProj(2, nR+4);
      PyQuICC::CoreWrapper::fillMatrix(tmpProj, pValue);
      Py_DECREF(pValue);

      // Finalize the Python wrapper
      PyQuICC::CoreWrapper::finalize();

      // Store integral projector
      SparseMatrix temp = tmpProj*tmpI1*tmpR1*tmpR2;
      this->mIntgOp  = (temp.row(0)-temp.row(1));
      IVariableAsciiWriter::init();
   }

   void ShellTorPolUniformVorticityWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
      assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

      MHDFloat ro = this->mPhysical.find(NonDimensional::Ro::id())->second->value();
      MHDFloat ri = ro*this->mPhysical.find(NonDimensional::RRatio::id())->second->value();

      // Compute integral over Chebyshev expansion and sum harmonics
      this->mUVx = 0.0;
      this->mUVy = 0.0;
      this->mUVz = 0.0;

      #ifdef QUICC_SPATIALSCHEME_SLFM
         double factor = 1.0;
         // Loop over harmonic order m
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
        	int m = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			if(m > 1){
				continue;
			}

            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

               if( l!=1){
            	   continue;
               }

				if(m==0){
					this->mUVz = this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real());
				} else {
					this->mUVx = - this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real())*std::sqrt(2);
					this->mUVy = this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).imag())*std::sqrt(2);
				}
            }
         }
      #endif //defined QUICC_SPATIALSCHEME_SLFM
      #ifdef QUICC_SPATIALSCHEME_SLFL
         // Loop over harmonic degree l
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if( l!=1){
         	   continue;
            }

            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){
            	int m = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);


            	if(m==0){
            		this->mUVz = this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real());
            	} else {
            		this->mUVx = - this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real())*std::sqrt(2);
            		this->mUVy = this->mIntgOp.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).imag())*std::sqrt(2);
            	}
            }
         }
      #endif //QUICC_SPATIALSCHEME_SLFL

      // Normalize the integral with int_ri^ror^4dr minus the boundaries
      this->mUVz /= this->mVolume;
      this->mUVx /= this->mVolume;
      this->mUVy /= this->mVolume;
   }

   void ShellTorPolUniformVorticityWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array vorticity(3);

         vorticity(0) = this->mUVx;
         vorticity(1) = this->mUVy;
         vorticity(2) = this->mUVz;

         MPI_Allreduce(MPI_IN_PLACE, vorticity.data(), 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mUVx = vorticity(0);
         this->mUVy = vorticity(1);
         this->mUVz = vorticity(2);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(QuICCEnv().allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mUVx << '\t' << this->mUVy << '\t' << this->mUVz << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mUVz) || std::isnan(this->mUVx) || std::isnan(this->mUVy))
      {
         QuICCEnv().abort("Uniform Vorticity is NaN!");
      }
   }

}
}
}
