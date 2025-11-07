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

protected:
private:
};

} // namespace DenseSM
} // namespace QuICC


#endif // define QUICC_DENSESM_IGENERICPROFILE_HPP
