#ifndef included_NodalVariable
#define included_NodalVariable

#include "vectors/Variable.h"


namespace AMP {
namespace Discretization {


/**
 * \class NodalVariable
 * \brief A class used to describe a simple nodal variable
 */
class NodalVariable : public AMP::LinearAlgebra::Variable
{
public:
    //! Default constructor for a nodal variable
    NodalVariable( int DOFsPerNode, const std::string &name );

    //!  Return the number of DOFs required per node
    virtual size_t 	DOFsPerObject () const;

    //! Return the variable ID type (not an ID of the variable)
    virtual  size_t  variableID() const;

private:
    
    // Number of DOFs per node
    size_t d_DOFsPerNode;

};


}
}

#endif


