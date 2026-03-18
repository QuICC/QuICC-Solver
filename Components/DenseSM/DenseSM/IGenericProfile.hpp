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

   /**
    * @brief Evaluate function on grid, non-multiprecision version
    */
   Array evaluateLP(const Array& r, const int l, const int m) const
   {
      // Cast input r (single-precision) to Internal::Array (multiprecision)
      Internal::Array rMP = r.template cast<Internal::MHDFloat>().eval();

      // Evaluate in multiprecision, then cast result back to single-precision Array
      return this->evaluate(rMP, l, m).template cast<MHDFloat>().eval();
   }

   /**
     * @brief Get the name of the profile
     */
    virtual const std::string& getName() const { return mName; }

protected:

   /**
     * @brief Name of the profile (for output)
     */
    std::string mName = "";

private:
};

} // namespace DenseSM
} // namespace QuICC


#endif // QUICC_DENSESM_IGENERICPROFILE_HPP
