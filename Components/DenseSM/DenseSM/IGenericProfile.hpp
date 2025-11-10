/**
 * @file IGenericProfile.hpp
 * @brief Implementation of the generic interface for a profile function (e.g. radial profiles)
 */

#ifndef QUICC_DENSESM_IGENERICPROFILE_HPP
#define QUICC_DENSESM_IGENERICPROFILE_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Typedefs.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

/**
 * @brief Implementation of the generic profile function
 */
class IGenericProfile
{
public:
   /**
    * @brief Constructor
    */
   IGenericProfile() = default;

   /**
    * @brief Destructor
    */
   virtual ~IGenericProfile() = default;

   /**
    * @brief Required spectral truncation
    */
   virtual int nN() const = 0;

   /**
    * @brief Nonzero harmonic degrees
    */
   virtual std::vector<int> ls() const = 0;

   /**
    * @brief Evaluate function on grid
    */
   virtual Internal::Array evaluate(const Internal::Array& r, const int l,
      const int m) const = 0;

protected:


protected:
private:
};

} // namespace DenseSM
} // namespace QuICC


#endif // QUICC_DENSESM_IGENERICPROFILE_HPP
