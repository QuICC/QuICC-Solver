/**
 * @file ModelChecker.h Check compatiblity between model and framework
 */

#ifndef QUICC_VERSION_FRAMEWORK_MODELCHECKER_H
#define QUICC_VERSION_FRAMEWORK_MODELCHECKER_H

#include "QuICC/Version/Framework.hpp"

// Check model version information is defined
//
#ifndef QUICC_VERSION_MODEL_MAJOR
   #error "QUICC_VERSION_MODEL_MAJOR is not defined!"
#endif

#ifndef QUICC_VERSION_MODEL_MINOR
   #error "QUICC_VERSION_MODEL_MINOR is not defined!"
#endif

#ifndef QUICC_VERSION_MODEL_PATCH
   #error "QUICC_VERSION_MODEL_PATCH is not defined!"
#endif

// Check MAJOR version compatiblity
//
#if QUICC_VERSION_FRAMEWORK_MAJOR != QUICC_VERSION_MODEL_MAJOR
   #error "Model and QuICC major version number don't match!"
#endif //QUICC_VERSION_FRAMEWORK_MAJOR != QUICC_VERSION_MODEL_MAJOR

// Check MINOR version compatiblity
//
#if QUICC_VERSION_FRAMEWORK_MINOR > QUICC_VERSION_MODEL_MINOR
   #warning "QuICC minor version is more recent than Model!"
#elif QUICC_VERSION_FRAMEWORK_MINOR < QUICC_VERSION_MODEL_MINOR
   #error "Model minor version is more recent than QuICC!"
#endif //QUICC_VERSION_FRAMEWORK_MINOR > QUICC_VERSION_MODEL_MINOR

#endif // QUICC_VERSION_FRAMEWORK_MODELCHECKER_H
