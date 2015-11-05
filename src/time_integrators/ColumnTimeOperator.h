#ifndef included_AMP_ColumnTimeOperator
#define included_AMP_ColumnTimeOperator

#include "operators/ColumnOperator.h"

namespace AMP{
namespace TimeIntegrator{

  /**
   * The ColumnTimeOperator class is meant to handle the case where
   * the TimeOperator is a multi-physics operator with ColumnOperator's
   * representing the mass and rhs operators. It is meant to allow for
   * the use of ColumnSolvers in block diagonal or block triangular preconditioners
   * for Newton-Krylov methods in implicit time integration. This operator is currently
   * derived from the ColumnOperator and not from the TimeOperator class. This results
   * in some duplication/reproduction of interfaces that would normally belong in the
   * TimeOperator class. The alternative is to allow multiple inheritance which is left
   * open as a possible design modification at this point.

   @see ColumnOperator
   @see TimeOperator
   */
  
class  ColumnTimeOperator: public AMP::Operator::ColumnOperator
  {
  public:

    /**
     * Main constructor for users. Expects a TimeOperatorParameters object.
     1. name: bLinearMassOperator
     description: boolean to indicate whether the mass operator is a linear operator, currently
     either both mass and rhs operators have to be linear or both have to be nonlinear
     type: bool
     default value: FALSE
     valid values: (TRUE, FALSE)
     optional field: yes
     
     2. name: bLinearRhsOperator
     description: boolean to indicate whether the rhs operator is a linear operator, currently
     either both mass and rhs operators have to be linear or both have to be nonlinear
     type: bool
     default value: FALSE
     valid values: (TRUE, FALSE)
     optional field: yes
     
     3. name: algebraicComponent
     description: for a DAE system this ColumnTimeOperator would have algebraic components. This field
     indicates the component that is an algebraic component, ie, having no time derivative. Note the limitation
     (temporary) to one algebraic component.
     type: integer
     default value: -1
     valid values: integer values
     optional field: yes, if no algebraic components are present
     
     4. name: CurrentDt
     description: current time step value
     type: double
     default value: 1.0e-08
     valid values: (positive real values)
     optional field: no

     5. name: CurrentTime
     description: current time
     type: double
     default value: 0.0
     valid values: ( real values)
     optional field: yes

     6. name: ScalingFactor
     description: 
     type: double
     default value:
     valid values:
     optional field: 
     
    */

    ColumnTimeOperator(AMP::shared_ptr<AMP::Operator::OperatorParameters > params);

    /**
     * virtual destructor
     */
    virtual ~ColumnTimeOperator();

    /**
     * This function is useful for re-initializing an operator
     * \param params
     *        parameter object containing parameters to change
     */
    virtual void reset(const AMP::shared_ptr<AMP::Operator::OperatorParameters>& params);
    
    /**
     * The apply routine for the column operator calls apply on each of the component operators
     */
    void apply(AMP::LinearAlgebra::Vector::const_shared_ptr u, 
	       AMP::LinearAlgebra::Vector::shared_ptr f ) override;
    
    /**
     * \param op
     *            shared pointer to an operator to append to the existing column of operators
     */
    virtual void append(AMP::shared_ptr< AMP::Operator::Operator > op);
    
    /**
     * This function registers a rhs operator with the TimeOperator class
     @param [in] op : shared pointer to Operator
    */
    void registerRhsOperator(AMP::shared_ptr< AMP::Operator::ColumnOperator > op) {d_pRhsOperator = op; }

    /**
     * This function registers a mass operator with the TimeOperator class. Not necessary
     * for FD or FVM discretizations
     */
   void registerMassOperator(AMP::shared_ptr< AMP::Operator::ColumnOperator > op) {d_pMassOperator = op; }
    
   /**
    * return a shared pointer to the rhs operator
    */
    AMP::shared_ptr< Operator > getRhsOperator(void){ return d_pRhsOperator; }
    
    /**
     * return a shared pointer to the mass operator
     */
    AMP::shared_ptr< Operator > getMassOperator(void){ return d_pMassOperator; }
    
  protected:
    
    ColumnTimeOperator();
    
    void getFromInput(const AMP::shared_ptr<AMP::Database>&);
    
    /**
     * rhs and mass operators are intentionally chosen to be column
     */
    AMP::shared_ptr< AMP::Operator::ColumnOperator > d_pRhsOperator;
    
    AMP::shared_ptr< AMP::Operator::ColumnOperator > d_pMassOperator;
    
    AMP::shared_ptr<AMP::LinearAlgebra::Vector>  d_pPreviousTimeSolution;
    
    AMP::shared_ptr<AMP::LinearAlgebra::Vector>  d_pSourceTerm;
    
  private:
    
  bool d_bCreateLinearTimeOperators;

  int d_iAlgebraicComponent;
  
  double d_dCurrentDt;
};

}
}

#endif
