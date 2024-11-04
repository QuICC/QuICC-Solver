/**
 * @file IPlaneModelBackend.hpp
 * @brief Base model backend for plane models
 */

#ifndef QUICC_MODEL_IPLANEMODELBACKEND_HPP
#define QUICC_MODEL_IPLANEMODELBACKEND_HPP

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
 * @brief Base model backend for Plane model
 */
class IPlaneModelBackend : public ICppModelBackend
{
public:
   /**
    * @brief Constructor
    */
   IPlaneModelBackend() = default;

   /**
    * @brief Destructor
    */
   virtual ~IPlaneModelBackend() = default;

protected:
private:
};

} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_IPLANEMODELBACKEND_HPP
