/**
 * @file ShellDipolarityWriter.hpp
 * @brief Implementation of the dipolarity in a spherical shell
 */

#ifndef QUICC_IO_VARIABLE_SHELLDIPOLARITYWRITER_HPP
#define QUICC_IO_VARIABLE_SHELLDIPOLARITYWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the dipolarity in a shell
 */
class ShellDipolarityWriter : public IVariableAsciiWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   ShellDipolarityWriter(const std::string& prefix, const std::string& type);

   /**
    * @brief Destructor
    */
   ~ShellDipolarityWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

   /**
    * @brief Compute dipolarity of magnetic field
    */
   void compute(Transform::TransformCoordinatorType& coord);

   /**
    * @brief Requires heavy calculation?
    */
   virtual bool isHeavy() const;

   /**
    * @brief Set truncation of CMB spectrum output
    */
   void setCmbTruncation(const int nL);

protected:
   /**
    * @brief Prepare spectral field data for computation
    */
   void prepareInput(const FieldComponents::Spectral::Id sId,
      Transform::TransformCoordinatorType& coord);

   /**
    * @brief Write State to file
    */
   virtual void writeContent();

   /**
    * @brief Data ordering is m slowest
    */
   bool mHasMOrdering;

private:
   /**
    * @brief Reset energy storage
    */
   void resetEnergy();

   /**
    * @brief CMB spectrum output truncation
    */
   int mCmbNl;

   /**
    * @brief Dipolarity
    */
   MHDFloat mDipolarity;

   /**
    * @brief Axial dipole component
    */
   MHDFloat mAxialDipole;

   /**
    * @brief Non axial dipole component
    */
   MHDComplex mNonAxialDipole;

   /**
    * @brief Energy spectrum at CMB
    */
   Array mCmbSpectrum;

   /**
    * @brief Boundary value operators
    */
   std::map<int, Array> mValue;
};

/// Typedef for a shared pointer to and DipolarityWriter
typedef std::shared_ptr<ShellDipolarityWriter> SharedShellDipolarityWriter;

inline bool ShellDipolarityWriter::isHeavy() const
{
   return true;
}

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SHELLDIPOLARITYWRITER_HPP
