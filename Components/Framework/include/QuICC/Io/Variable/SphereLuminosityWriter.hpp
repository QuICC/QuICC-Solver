/**
 * @file SphereLuminosityWriter.hpp
 * @brief Implementation of the Luminosity in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERELUMINOSITYWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERELUMINOSITYWRITER_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Io/Variable/IVariableAsciiWriter.hpp"
#include "DenseSM/Worland/RadialTorPolFunction.hpp"


namespace QuICC {

namespace Io {

namespace Variable {

/**
 * @brief Implementation of the Nusselt number in a sphere
 */
class SphereLuminosityWriter : public IVariableAsciiWriter
{
public:
   /**
    * @brief Constructor
    *
    * @param prefix Prefix to use for file name
    * @param type Type of the file (typically scheme name)
    */
   SphereLuminosityWriter(const std::string& prefix, const std::string& type, std::vector<std::shared_ptr<QuICC::DenseSM::Worland::RadialTorPolFunction>> pF);

   /**
    * @brief Destructor
    */
   virtual ~SphereLuminosityWriter() = default;

   /**
    * @brief Initialise the operator, transform and file
    */
   virtual void init();

   /**
    * @brief Requires heavy calculation?
    */
   virtual bool isHeavy() const;

protected:
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
    * @brief Luminosity
    */
   MHDFloat mLuminosity;

   /**
    * @brief Nusselt number
    */
   MHDFloat mNusselt;

   /*
    * @brief Spherical volume to normalize energy to energy density
    */
   MHDFloat mSb;

   /**
    * @brief Origin projector
    */
   Matrix mBoundary;

   /**
    * @brief shared pointers to density*Temperature*kappa profile 
    */
   std::shared_ptr<QuICC::DenseSM::Worland::RadialTorPolFunction> mpRhoTempKappa;

   /**
    * @brief shared pointers to D1ConductiveEntropy profile 
    */
   std::shared_ptr<QuICC::DenseSM::Worland::RadialTorPolFunction> mpD1Sc;
};

/// Typedef for a shared pointer of a HDF5 state file writer
typedef std::shared_ptr<SphereLuminosityWriter> SharedSphereLuminosityWriter;

inline bool SphereLuminosityWriter::isHeavy() const
{
   return false;
}

} // namespace Variable
} // namespace Io
} // namespace QuICC

#endif // QUICC_IO_VARIABLE_SPHERELUMINOSITYWRITER_HPP
