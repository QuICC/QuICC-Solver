/**
 * @file ShellTorPolTorqueWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a spherical shell
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
#include "QuICC/Io/Variable/ShellTorPolTorqueWriter.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/Ro.hpp"
#include "QuICC/NonDimensional/RRatio.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Io/Variable/EnergyTags.hpp"
#include "QuICC/TypeSelectors/ScalarSelector.hpp"
#include "QuICC/PyQuICC/CoreWrapper.hpp"
#include<typeinfo>
#include<iostream>
#include<numpy/npy_3kcompat.h>
#include<math.h>

namespace QuICC {

namespace Io {

namespace Variable {

   ShellTorPolTorqueWriter::ShellTorPolTorqueWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL, IVariableAsciiWriter::EXTEND), mTorque(1.2345),mProj(),mFactor(0.)
   {

   }

   ShellTorPolTorqueWriter::~ShellTorPolTorqueWriter()
   {
   }

   void ShellTorPolTorqueWriter::init()
   {

	   // initialize for every core the compute flag to false
	   this->mComputeFlag = false;
	   #ifdef QUICC_SPATIALSCHEME_SLFM
	   // Loop over harmonic order m
	   for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
	   {
		   // determine current m
		   int m = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
		   // m = 0, no factor of two


		   for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
		   {
			   int l = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

			   if( m==0 && l==1){
				   this->mComputeFlag = true;
				   break;
			   }
		   }
	   }

       #endif //defined QUICC_SPATIALSCHEME_SLFM
	   #ifdef QUICC_SPATIALSCHEME_SLFL
	   // Loop over harmonic degree l
	   for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
	   {
		   int l = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);


		   for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){
			   int m = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

			   if( m ==0 && l==1){
				   this->mComputeFlag = true;
				   break;
			   }
		   }

	   }
       #endif //QUICC_SPATIALSCHEME_SLFL

	   // if necessary prepare projection operators and weights
	   if(this->mComputeFlag){

		   // Spherical shell volume: 4/3*pi*(r_o^3 - r_i^3)
		   MHDFloat ro = this->mPhysical.find(NonDimensional::Ro::id())->second->value();
		   MHDFloat ri = ro*this->mPhysical.find(NonDimensional::RRatio::id())->second->value();
		   MHDFloat E =this->mPhysical.find(NonDimensional::Ekman::id())->second->value();

		   double pi = 3.14159265358979;
		   this->mFactor = 8./3.*pi*ri*ri*E;

		   // Initialise python wrapper
		   PyQuICC::CoreWrapper::init();
		   PyQuICC::CoreWrapper::import("quicc.geometry.spherical.shell_radius");

		   // Prepare arguments
		   PyObject *pArgs, *pValue;
		   pArgs = PyTuple_New(4);

		   // ... compute a, b factors
		   PyObject *pTmp = PyTuple_New(2);
		   PyTuple_SetItem(pTmp, 0, PyFloat_FromDouble(ro));
		   PyTuple_SetItem(pTmp, 1, PyFloat_FromDouble(this->mPhysical.find(NonDimensional::RRatio::id())->second->value()));
		   PyQuICC::CoreWrapper::setFunction("linear_r2x");
		   pValue = PyQuICC::CoreWrapper::callFunction(pTmp);
		   MHDFloat a =PyFloat_AsDouble(PyTuple_GetItem(pValue,0));
		   MHDFloat b =PyFloat_AsDouble(PyTuple_GetItem(pValue,1));

		   // prepare the next function call
		   pArgs = PyTuple_New(2);

		   // create boundray conditions (none)
		   PyObject * pDict1;
		   pDict1 = PyDict_New();

		   PyDict_SetItem(pDict1, PyLong_FromLong(0), PyLong_FromLong(20));
		   PyTuple_SetItem(pArgs, 1, pDict1);

		   // Get resolution
		   //int nR = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>() + 2;
		   int nR = this->res().sim().dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

		   PyTuple_SetItem(pArgs, 0, PyLong_FromLong(nR));

		   // Call zblk
		   PyQuICC::CoreWrapper::setFunction("zblk");
		   pValue = PyQuICC::CoreWrapper::callFunction(pArgs);
		   // Fill matrix and cleanup
		   SparseMatrix valueProj(nR,nR);
		   PyQuICC::CoreWrapper::fillMatrix(valueProj, pValue);
		   Py_DECREF(pValue);

		   // add specifics to boundary conditions
		   pDict1 = PyDict_New();

		   PyObject* pDict2;
		   pDict2 = PyDict_New();
		   PyDict_SetItem(pDict2, PyUnicode_FromString("a"), PyFloat_FromDouble(a));
		   PyDict_SetItem(pDict2, PyUnicode_FromString("b"), PyFloat_FromDouble(b));

		   // insert pDict2 into pDict1
		   PyDict_SetItem(pDict1, PyLong_FromLong(0), PyLong_FromLong(21));
		   PyDict_SetItem(pDict1, PyUnicode_FromString("c"), pDict2);
		   PyTuple_SetItem(pArgs, 1, pDict1);
		   PyTuple_SetItem(pArgs, 0, PyLong_FromLong(nR));

		   // Call avg
		   PyQuICC::CoreWrapper::setFunction("zblk");
		   pValue = PyQuICC::CoreWrapper::callFunction(pArgs);

		   // Fill matrix and cleanup
		   SparseMatrix diffProj(nR,nR);
		   PyQuICC::CoreWrapper::fillMatrix(diffProj, pValue);
		   Py_DECREF(pValue);

		   // select only the correct boundary condition
		   this->mProj = (diffProj*ri - valueProj).row(1);
		   PyQuICC::CoreWrapper::finalize();
	   }

	   // call init in the base class
	   IVariableAsciiWriter::init();
   }

   void ShellTorPolTorqueWriter::compute(Transform::TransformCoordinatorType& coord)
   {
	   this->mTorque = 0;
	   if(this->mComputeFlag){

		   // get iterator to field
		   vector_iterator vIt;
		   vector_iterator_range vRange = this->vectorRange();
		   assert(std::distance(vRange.first, vRange.second) == 1);
		   assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
		   assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

		   // retrieve ri
		   MHDFloat ro = this->mPhysical.find(NonDimensional::Ro::id())->second->value();
		   MHDFloat ri = ro*this->mPhysical.find(NonDimensional::RRatio::id())->second->value();



#ifdef QUICC_SPATIALSCHEME_SLFM
		   // Loop over harmonic order m
		   for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
		   {
			   // determine current m
			   int m = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			   for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
			   {
				   int l = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				   if(l==1 && m==0){

					   // compute torque
					   MHDFloat temp = this->mProj.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real());
					   this->mTorque = temp*this->mFactor;
					   break;
				   }

			   }
		   }
#endif //defined QUICC_SPATIALSCHEME_SLFM

#ifdef QUICC_SPATIALSCHEME_SLFL
		   // Loop over harmonic degree l
		   for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
		   {
			   // determine current l
			   int l = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);

			   // loop over order m
			   for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){

				   // determine current m
				   int m = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);

				   if(l==1 && m==0){

					   // compute torque
					   MHDFloat temp = this->mProj.dot(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR).slice(k).col(j).real());
					   this->mTorque = temp*this->mFactor;
					   break;
				   }
			   }

		   }
#endif //QUICC_SPATIALSCHEME_SLFL

	   }

   }

   void ShellTorPolTorqueWriter::write()
   {
	   // Create file
	   this->preWrite();


	   // message pass the torque
	   #ifdef QUICC_MPI

	   MHDFloat Torque = this->mTorque;

	   MPI_Allreduce(MPI_IN_PLACE, &Torque, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	   this->mTorque = Torque;
	   #endif //QUICC_MPI


	   // Check if the workflow allows IO to be performed
	   if(QuICCEnv().allowsIO())
	   {
		   this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTorque << std::endl;
	   }

	   // Close file
	   this->postWrite();

	   // Abort if kinetic energy is NaN
	   if(std::isnan(this->mTorque))
	   {
		   QuICCEnv().abort("Toroidal torque is NaN!");
	   }
   }

}
}
}
