#ifndef included_AMP_ScalarN2GZAxisMap
#define included_AMP_ScalarN2GZAxisMap

#include "operators/map/Map3to1to3.h"
#include "operators/map/Map3to1to3Parameters.h"
#include "discretization/createLibmeshElements.h"

/* Libmesh files */
#include "fe_type.h"
#include "fe_base.h"
#include "elem.h"
#include "quadrature.h"

namespace AMP {
namespace Operator {


typedef AMP::Operator::Map3to1to3Parameters  ScalarN2GZAxisMapParameters;


/**
 * \class  ScalarN2GZAxisMap
 * \brief  A class used to reduce a 3D problem to 1D, transfer the solution, and map back to 3D
 * \details  This class inherites from Map3to1to3, and performs a reduction from 3D to 1D, 
 *    transfers the solution, then maps back to 3D.  It accomplishes this by taking the average
 *    value for each point along the z-axis, creating a 1D function of z, then transfering that 
 *    solution.  To behave correctly, the different nodes must be aligned in the z-direction.
 */
class ScalarN2GZAxisMap : public AMP::Operator::Map3to1to3
{
public:

    /** \brief  Returns true if MapType = "ScalarN2GZAxis"
     * \param[in] s  A string extracted from the MapType line in a MeshToMeshMap db
     * \return  True iff s == "ScalarN2GZAxis"
     */
    static bool validMapType ( const std::string &s );

    /** \brief  Typedef to identify the parameters class of this operator
     */
    typedef  ScalarN2GZAxisMapParameters   Parameters;

    //!  The base tag used in communication.
    enum { CommTagBase = 20000 };

    /** \brief   Standard constructor
     * \param[in] params  Input parameters
     */
    ScalarN2GZAxisMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &params );

    //! Destructor
    virtual ~ScalarN2GZAxisMap ();

protected:
    // Implimented buildMap routine
    virtual std::multimap<double,double>  buildMap( const AMP::LinearAlgebra::Vector::shared_ptr, 
        const AMP::Mesh::MeshIterator& );

    // Implimented buildReturn routine
    virtual void buildReturn( AMP::LinearAlgebra::Vector::shared_ptr, 
        const AMP::Mesh::MeshIterator&, const std::multimap<double,double>& );

private:

    boost::shared_ptr < ::FEBase > d_fe; 
    Discretization::createLibmeshElements libmeshElements;

    void createSurfaceFEBase();
};


} // Operator namespace
} // AMP namespace

#endif
