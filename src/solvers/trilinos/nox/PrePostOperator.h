#ifndef included_AMP_ThermalPrePostOperator
#define included_AMP_ThermalPrePostOperator

#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Solver_Generic.H"

#include "utils/shared_ptr.h"
#include "vectors/Vector.h"


namespace AMP {
namespace Solver {


/**
 * \class ThermalPrePostOperator
 * \brief A class used to interact with the NOX solver for the MulitPelletThermalSolver
 */
class PrePostOperator : public NOX::Abstract::PrePostOperator
{
public:
    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a mesh manager object.
     */
    typedef AMP::shared_ptr<AMP::Solver::PrePostOperator> shared_ptr;

    // Destructor
    virtual ~PrePostOperator() {}

    // Functions inherited from NOX::Abstract::PrePostOperator
    virtual void runPreIterate( const NOX::Solver::Generic &solver )  = 0;
    virtual void runPostIterate( const NOX::Solver::Generic &solver ) = 0;
    virtual void runPreSolve( const NOX::Solver::Generic &solver )    = 0;
    virtual void runPostSolve( const NOX::Solver::Generic &solver )   = 0;

    // Functions for pre/post apply
    virtual void runPreApply( AMP::LinearAlgebra::Vector::const_shared_ptr x,
                              AMP::LinearAlgebra::Vector::shared_ptr f,
                              bool exact ) = 0;
    virtual void runPostApply( AMP::LinearAlgebra::Vector::const_shared_ptr x,
                               AMP::LinearAlgebra::Vector::shared_ptr f,
                               bool exact ) = 0;

protected:
    // Constructor
    PrePostOperator() {}

private:
};
}
}

#endif
