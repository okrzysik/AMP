#ifndef included_AMP_Map3Dto1D
#define included_AMP_Map3Dto1D


#include "discretization/createLibmeshElements.h"
#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "operators/map/MapOperator.h"
#include "operators/map/MapOperatorParameters.h"
#include "utils/shared_ptr.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"

#include <string>

namespace AMP {
namespace Operator {

/**
 *  Class Map3Dto1D is used to map the solution on the boundary of the 3D mesh to
 *  the 1D location by interpolation. The 3D mesh is the d_MeshAdapter from the
 *  Operator class.
 */
class Map3Dto1D : public MapOperator {
public:
    //!  Constructor
    explicit Map3Dto1D( const AMP::shared_ptr<OperatorParameters> &params );

    //!  Destructor
    virtual ~Map3Dto1D() {}

    /**
      This function reads the entries of the database for the operator
      and can also be used to change the parameters if required.
     */
    void reset( const AMP::shared_ptr<OperatorParameters> & );

    /**
      For this operator the apply function would map the solution by interpolation from
      NodalScalar Vector in u to Simple Vector in r.
      @param [in]  u input vector.
      @param [out] r output vector.
     */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    void apply_Gauss( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                      AMP::LinearAlgebra::Vector::shared_ptr f );

    void apply_Nodal( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                      AMP::LinearAlgebra::Vector::shared_ptr f );

    AMP::LinearAlgebra::Variable::shared_ptr createInputVariable( const std::string &name,
                                                                  int = -1 )
    {
        return d_inpVariable->cloneVariable( name );
    }

    AMP::LinearAlgebra::Variable::shared_ptr createOutputVariable( const std::string &name,
                                                                   int = -1 )
    {
        return d_outVariable->cloneVariable( name );
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() { return d_inpVariable; }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() { return d_outVariable; }

    /**
      This member function is used to set the 1D Z locations std vector.
      */
    void setZLocations( std::vector<double> zloc ) { d_zLocations = zloc; }

    void makeZLocationsConsistent();

    void setVector( AMP::LinearAlgebra::Vector::shared_ptr vec )
    {
        outputVec = vec->subsetVectorForVariable( d_outVariable );
    }

protected:
    AMP::LinearAlgebra::Vector::shared_ptr outputVec;

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_inpVariable;

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_outVariable;

    std::vector<double> d_zLocations; /**< std vector to store 1D z locations. */

    bool d_useGaussVec;

private:
    Discretization::createLibmeshElements libmeshElements;
};
}
}

#endif
