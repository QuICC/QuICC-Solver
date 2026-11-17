/**
 * @file TestBackend.hpp
 * @brief Test model backend
 */

#ifndef QUICC_TESTSUITE_FRAMEWORK_IO_TESTBACKEND_HPP
#define QUICC_TESTSUITE_FRAMEWORK_IO_TESTBACKEND_HPP

// System includes
//
#include <map>
#include <memory>
#include <string>
#include <vector>

// Project includes
//
#include "QuICC/Model/ISphericalModelBackend.hpp"

namespace QuICC {

namespace TestSuite {

namespace Framework {

namespace Io {

/**
 * @brief Interface for model backend
 */
class TestBackend : public QuICC::Model::ISphericalModelBackend
{
public:
   /**
    * @brief Constructor
    */
   TestBackend(const  std::string& scheme);

   /**
    * @brief Destructor
    */
   virtual ~TestBackend() = default;

   /**
    * @brief Get vector of names for the physical fields
    */
   virtual std::vector<std::string> fieldNames() const override;

   /**
    * @brief Get vector of names for the nondimensional parameters
    */
   virtual std::vector<std::string> paramNames() const override;

   /**
    * @brief Get vector of bools about periodic box
    */
   virtual std::vector<bool> isPeriodicBox() const override;

   /**
    * @brief Get automatically computed parameters based on input parameters
    *
    * @param cfg  Input parameters
    */
   virtual std::map<std::string, MHDFloat> automaticParameters(
      const std::map<std::string, MHDFloat>& cfg) const override;

   /**
    * @brief Get equation information
    *
    * @param info Equation information
    * @param fId  Field ID
    * @param res  Resolution object
    */
   virtual void equationInfo(QuICC::Model::EquationInfo& info, const SpectralFieldId& fId,
      const Resolution& res) const override;

   /**
    * @brief Get operator information
    *
    * @param info       Equation information
    * @param fId        Field ID
    * @param res        Resolution object
    * @param coupling   Equation/Field coupling information
    * @param bcs        Boundary conditions
    */
   virtual void operatorInfo(QuICC::Model::OperatorInfo& info, const SpectralFieldId& fId,
      const Resolution& res, const Equations::Tools::ICoupling& coupling,
      const BcMap& bcs, const bool allowGalerkin) const override;

   /**
    * @brief Build model matrix
    *
    * @param rModelMatrix  Input/Output matrix to fill with operators
    * @param opId          Type of model matrix
    * @param imRange       Coupled fields
    * @param matIdx        Matrix index
    * @param bcType        Boundary condition scheme (Tau vs Galerkin)
    * @param res           Resolution object
    * @param eigs          Indexes of other dimensions
    * @param bcs           Boundary conditions
    * @param nds           Nondimensional parameters
    */
   virtual void modelMatrix(DecoupledZSparse& rModelMatrix,
      const std::size_t opId,
      const Equations::CouplingInformation::FieldId_range imRange,
      const int matIdx, const std::size_t bcType, const Resolution& res,
      const std::vector<MHDFloat>& eigs, const BcMap& bcs,
      const NonDimensional::NdMap& nds) const override;

   /**
    * @brief Build galerkin stencil
    *
    * @param mat     Input/Output matrix to fill with stencil
    * @param fId     Field ID
    * @param matIdx  Matrix index
    * @param res     Resolution object
    * @param eigs          Indexes of other dimensions
    * @param makeSquare Truncate stencil to obtain square matrix?
    * @param bcs           Boundary conditions
    * @param nds           Nondimensional parameters
    */
   virtual void galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fId,
      const int matIdx, const Resolution& res,
      const std::vector<MHDFloat>& eigs, const bool makeSquare,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

   /**
    * @brief Build explicit block
    *
    * @param rModelMatrix  Input/Output matrix to fill with operators
    * @param fId           Equation field ID
    * @param opId          Type of explicit operator
    * @param fieldId       Coupled field ID
    * @param matIdx        Matrix index
    * @param res           Resolution object
    * @param eigs          Indexes of other dimensions
    * @param bcs           Boundary conditions
    * @param nds           Nondimensional parameters
    */
   virtual void explicitBlock(DecoupledZSparse& mat, const SpectralFieldId& fId,
      const std::size_t opId, const SpectralFieldId fieldId, const int matIdx,
      const Resolution& res, const std::vector<MHDFloat>& eigs,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const override;

protected:
   /**
    * @brief Number of boundary conditions
    *
    * @param fId  Field ID
    */
   int nBc(const SpectralFieldId& fId) const override;

   /**
    * @brief Apply tau line for boundary condition
    *
    * @param mat     Input/Output matrix to apply tau line to
    * @param rowId   ID of field of equation
    * @param colId   ID of field
    * @param l       Harmonic degree
    * @param opts    Options
    * @param nN      1D dimension
    * @param bcs     Boundary conditions
    * @param nds     Nondimensional parameters
    * @param isSplitOperator  Is second operator of split 4th order system?
    */
   void applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const int l,
      std::shared_ptr<QuICC::Model::details::BlockOptions> opts, const int nN,
      const BcMap& bcs, const NonDimensional::NdMap& nds,
      const bool isSplitOperator) const override {};

   /**
    * @brief Apply galerkin stencil for boundary condition
    *
    * @param mat     Input/Output matrix to apply stencil to
    * @param rowId   ID of field of equation
    * @param colId   ID of field
    * @param lr      Row space harmonic degree
    * @param lc      Column space harmonic degree
    * @param opts    Options
    * @param nN      1D dimension
    * @param bcs     Boundary conditions
    * @param nds     Nondimensional parameters
    */
   void applyGalerkinStencil(SparseMatrix& decMat, const SpectralFieldId& rowId,
      const SpectralFieldId& colId, const int lr, const int lc,
      std::shared_ptr<QuICC::Model::details::BlockOptions> opts, const int nNr, const int nNc,
      const BcMap& bcs, const NonDimensional::NdMap& nds) const override {};

private:
   /**
    * @brief Operators are complex?
    *
    * @param fId  Field ID
    */
   bool isComplex(const SpectralFieldId& fId) const final;

   /**
    * @brief Base 1D dimension
    *
    * @fId  Field ID
    */
   int baseNn(const int l, const Resolution& res) const;

   /**
    * @brief Get coupled fields
    *
    * @param fId  Field ID
    */
   SpectralFieldIds implicitFields(const SpectralFieldId& fId) const final;

   /**
    * @brief Scheme ID
    */
   std::string mScheme;

};

} // namespace Io
} // namespace Framework
} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_FRAMEWORK_IO_TESTBACKEND_HPP
