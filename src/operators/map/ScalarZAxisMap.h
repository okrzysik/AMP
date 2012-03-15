#ifndef included_AMP_ScalarZAxisMap
#define included_AMP_ScalarZAxisMap

#include "operators/map/Map3to1to3.h"
#include "operators/map/Map3to1to3Parameters.h"


namespace AMP {
namespace Operator {


typedef AMP::Operator::Map3to1to3Parameters  ScalarZAxisMapParameters;


/**
 * \class  ScalarZAxisMap
 * \brief  A class used to reduce a 3D problem to 1D, transfer the solution, and map back to 3D
 * \details  This class inherites from Map3to1to3, and performs a reduction from 3D to 1D, 
 *    transfers the solution, then maps back to 3D.  It accomplishes this by taking the average
 *    value for each point along the z-axis, creating a 1D function of z, then transfering that 
 *    solution.  To behave correctly, the different nodes must be aligned in the z-direction.
 */
class ScalarZAxisMap : public AMP::Operator::Map3to1to3
{
public:

    /** \brief  Returns true if MapType = "ScalarZAxis"
     * \param[in] s  A string extracted from the MapType line in a MeshToMeshMap db
     * \return  True iff s == "ScalarZAxis"
     */
    static bool validMapType ( const std::string &s );

    /** \brief  Typedef to identify the parameters class of this operator
     */
    typedef  ScalarZAxisMapParameters   Parameters;

    //!  The base tag used in communication.
    enum { CommTagBase = 20000 };

    /** \brief   Standard constructor
     * \param[in] params  Input parameters
     */
    ScalarZAxisMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> & );

    //! Destructor
    virtual ~ScalarZAxisMap ();

protected:
    // Implimented buildMap routine
    virtual std::multimap<double,double>  buildMap( const AMP::LinearAlgebra::Vector::shared_ptr, 
        const AMP::Mesh::Mesh::shared_ptr, const AMP::Mesh::MeshIterator& );

    // Implimented buildReturn routine
    virtual void buildReturn( AMP::LinearAlgebra::Vector::shared_ptr, const AMP::Mesh::Mesh::shared_ptr, 
        const AMP::Mesh::MeshIterator&, const std::multimap<double,double>& );
};


} // Operator namespace
} // AMP namespace

#endif
