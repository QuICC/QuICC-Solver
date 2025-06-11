/**
 * @file ISphericalModelBackend.hpp
 * @brief Base model backend for spherical models
 */

#ifndef QUICC_MODEL_ISPHERICALMODELBACKEND_HPP
#define QUICC_MODEL_ISPHERICALMODELBACKEND_HPP

// System includes
//
#include <map>
#include <memory>
#include <string>
#include <vector>

// Project includes
//
#include "QuICC/Model/ICppModelBackend.hpp"

namespace QuICC {

namespace Model {

/**
 * @brief Base model backend for RTC model
 */
class ISphericalModelBackend : public ICppModelBackend
{
public:
   /**
    * @brief Constructor
    */
   ISphericalModelBackend() = default;

   /**
    * @brief Destructor
    */
   virtual ~ISphericalModelBackend() = default;

protected:
private:
};

} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_ISPHERICALMODELBACKEND_HPP
