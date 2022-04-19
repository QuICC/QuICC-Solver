/**
 * @file CouplingInformation.hpp
 * @brief Implemenation of the coupling information of an equation
 */

#ifndef QUICC_EQUATIONS_COUPLINGINFORMATION_HPP
#define QUICC_EQUATIONS_COUPLINGINFORMATION_HPP

// Configuration includes
//

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"
#include "QuICC/Equations/Tools/ICoupling.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Implemenation of the coupling information of the equation
    */
   class CouplingInformation
   {
      public:
         /// Typedef to simplify notation for the field data
         typedef std::vector<SpectralFieldId> FieldIdVector;

         /// Typedef for an iterator for the field data
         typedef FieldIdVector::const_iterator  FieldId_iterator;

         /// Typedef for a range iterator for the field coupling data
         typedef std::pair<FieldId_iterator,FieldId_iterator>  FieldId_range;

         /**
          * @brief Enum for the equation type
          */
         enum EquationTypeId {
            /// Equation needs time marching
            PROGNOSTIC = 1,
            /// Equation needs a solver
            DIAGNOSTIC = 2,
            /// Equation is trivial (but with NL, QI, source, explicit)
            TRIVIAL = 3,
            /// Equation does nothing (ie "wrapper" for variable)
            WRAPPER = 4
         };

         /**
          * @brief Simple constructor
          */
         CouplingInformation();

         /**
          * @brief Simple empty destructor
          */
         virtual ~CouplingInformation();

         /**
          * @brief Get equation type
          */
         EquationTypeId equationType() const;

         /**
          * @brief Has a nonlinear term?
          */
         bool hasNonlinear() const;

         /**
          * @brief Has a quasi-inverse operator for nonlinear terms?
          */
         bool hasQuasiInverse() const;

         /**
          * @brief Has a source term?
          */
         bool hasSource() const;

         /**
          * @brief Has a boundary value?
          */
         bool hasBoundaryValue() const;

         /**
          * @brief Is the system complex?
          */
         bool isComplex() const;

         /**
          * @brief Is Galerkin system?
          */
         bool isGalerkin() const;

         /**
          * @brief Get index type
          */
         CouplingIndexType indexType() const;

         /**
          * @brief Get eigen tools
          */
         const Tools::ICoupling& couplingTools() const;

         /**
          * @brief Number of blocks in row of the system
          */
         int nBlocks() const;

         /**
          * @brief Number of systems
          */
         int nSystems() const;

         /**
          * @brief Size of Tau block
          *
          * @param idx System index
          */
         int tauN(const int idx) const;

         /**
          * @brief Size of Galerkin block
          *
          * @param idx System index
          */
         int galerkinN(const int idx) const;

         /**
          * @brief Index shift for the Galerkin block
          *
          * @param idx System index
          * @param dim dimension
          */
         int galerkinShift(const int idx, const int dim) const;

         /**
          * @brief Size of system
          *
          * @param idx System index
          */
         int systemN(const int idx) const;

         /**
          * @brief Index of the field within the system
          */
         int fieldIndex() const;

         /**
          * @brief Index of solver
          */
         int solverIndex() const;

         /**
          * @brief Start index in field values
          */
         int fieldStart() const;

         /**
          * @brief Number of RHS columns in solve
          *
          * @param idx System index
          */
         int rhsCols(const int idx) const;

         /**
          * @brief Add field to list of implicit timestep fields (coupled solve)
          *
          * @param fieldId Physical ID of the field
          * @param compId  Physical ID of the component
          */
         void addImplicitField(const std::size_t fieldId, const FieldComponents::Spectral::Id compId);

         /**
          * @brief Add field to list of explicit timestep fields
          *
          * @param fieldId Physical ID of the field
          * @param compId  Physical ID of the component
          */
         void addExplicitField(const std::size_t fieldId, const FieldComponents::Spectral::Id compId, const std::size_t opId);

         /**
          * @brief Sort the list of implicit fields and set field index accordingly
          *
          * @param fieldId Physical ID of the field
          * @param compId  Physical ID of the component
          */
         void sortImplicitFields(const std::size_t fieldId, const FieldComponents::Spectral::Id compId);

         /**
          * @brief Set settings for all systems
          *
          * @param typeId     Type of the solver
          * @param isComplex  Complex flag of solver
          * @param fieldStart Start index of the filed
          */
         void setGeneral(const CouplingInformation::EquationTypeId typeId, const bool isComplex, const int fieldStart);

         /**
          * @brief Set nonlinear flags
          *
          * @param hasNonlinear     Equation requires nonlinear computation?
          * @param hasQuasiInverse  Equation requires quasi-inverse on nonlinear terms?
          */
         void setNonlinear(const bool hasNonlinear, const bool hasQuasiInverse);

         /**
          * @brief Set source flag
          *
          * @param hasSource  Equation requires source term computation?
          */
         void setSource(const bool hasSource);

         /**
          * @brief Set boundary value flag
          *
          * @param hasBoundaryValue  Equation requires boundary value computation?
          */
         void setBoundaryValue(const bool hasBoundaryValue);

         /**
          * @brief Set system sizes
          *
          * @param nSystems         Number of systems
          * @param tauNs            Tau block size for each system
          * @param galerkinNs       Galerkin block size for each system
          * @param galerkinShifts   Galerkin index shifts
          * @param rhsCols          Number of columns in RHS for each system
          * @param systemNs         Number of columns in RHS for each system
          */
         void setSizes(const int nSystems, const ArrayI& tauNs, const ArrayI& galerkinNs, const MatrixI& galerkinShifts, const ArrayI& rhsCols, const ArrayI& systemNs);

         /**
          * @brief Set the index type
          *
          * @param id   Dimension type id of index
          * @param spCoupling   Coupling tools
          */
         void setIndexType(const CouplingIndexType id, Tools::SharedICoupling spCoupling);

         /**
          * @brief set the solver index
          */
         void setSolverIndex(const int idx);

         /**
          * @brief Get iterator to implicit fields
          */
         FieldId_range implicitRange() const;

         /**
          * @brief Get iterator to explicit fields
          */
         FieldId_range explicitRange(const std::size_t opId) const;

      protected:

      private:
         /**
          * @brief Storage for the implicit fields information
          */
         std::vector<SpectralFieldId>   mImplicitFields;

         /**
          * @brief Storage for the explicit linear fields information
          */
         std::vector<SpectralFieldId>   mExplicitLFields;

         /**
          * @brief Storage for the explicit nonlinear fields information
          */
         std::vector<SpectralFieldId>   mExplicitNLFields;

         /**
          * @brief Storage for the explicit nextstep fields information
          */
         std::vector<SpectralFieldId>   mExplicitNSFields;

         /**
          * @brief Storage for the equation type
          */
         EquationTypeId mEquationType;

         /**
          * @brief Storage for the nonlinear flag
          */
         bool mHasNonlinear;

         /**
          * @brief Storage for the quasi-inverse flag
          */
         bool mHasQuasiInverse;

         /**
          * @brief Storage for the source flag
          */
         bool mHasSource;

         /**
          * @brief Storage for the boundary value flag
          */
         bool mHasBoundaryValue;

         /**
          * @brief Storage for the complex flag
          */
         bool mIsComplex;

         /**
          * @brief Storage for the galerkin flag
          */
         bool mIsGalerkin;

         /**
          * @brief Storage for the index type
          */
         CouplingIndexType mIndexType;

         /**
          * @brief Storage for the number of systems
          */
         int mNSystems;

         /**
          * @brief Storage for the field index
          */
         int mFieldIndex;

         /**
          * @brief Storage for the solver index
          */
         int mSolverIndex;

         /**
          * @brief Storage for the field start index
          */
         int mFieldStart;

         /**
          * @brief Shared coupling tools
          */
         Tools::SharedICoupling mspCouplingTools;

         /**
          * @brief Storage for the Tau block sizes
          */
         ArrayI mTauNs;

         /**
          * @brief Storage for the Galerkin block sizes
          */
         ArrayI mGalerkinNs;

         /**
          * @brief Storage for the Galerkin index shifts
          */
         MatrixI mGalerkinShifts;

         /**
          * @brief Storage for the number of RHS in solver
          */
         ArrayI mRhsCols;

         /**
          * @brief Storage for the system sizes
          */
         ArrayI mSystemNs;
   };
}
}

#endif // QUICC_EQUATIONS_COUPLINGINFORMATION_HPP
