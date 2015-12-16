#ifndef included_AMP_ScalarN2GZAxisMap
#define included_AMP_ScalarN2GZAxisMap

#include "discretization/createLibmeshElements.h"
#include "operators/map/Map3to1to3.h"
#include "operators/map/Map3to1to3Parameters.h"


namespace AMP {
namespace Operator {


typedef AMP::Operator::Map3to1to3Parameters ScalarN2GZAxisMapParameters;


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
    static bool validMapType( const std::string &s );

    /** \brief  Typedef to identify the parameters class of this operator
     */
    typedef ScalarN2GZAxisMapParameters Parameters;

    //!  The base tag used in communication.
    enum { CommTagBase = 20000 };

    /** \brief   Standard constructor
     * \param[in] params  Input parameters
     */
    explicit ScalarN2GZAxisMap( const AMP::shared_ptr<AMP::Operator::OperatorParameters> &params );

    //! Destructor
    virtual ~ScalarN2GZAxisMap();

protected:
    // Implimented buildMap routine
    virtual std::multimap<double, double> buildMap( AMP::LinearAlgebra::Vector::const_shared_ptr,
                                                    const AMP::Mesh::Mesh::shared_ptr,
                                                    const AMP::Mesh::MeshIterator & );

    // Implimented buildReturn routine
    virtual void buildReturn( AMP::LinearAlgebra::Vector::shared_ptr,
                              const AMP::Mesh::Mesh::shared_ptr,
                              const AMP::Mesh::MeshIterator &,
                              const std::map<double, double> & );

    // Function to return the coordinates of the gauss points
    AMP::LinearAlgebra::Vector::const_shared_ptr getGaussPoints( const AMP::Mesh::MeshIterator & );

    // Internal vector with the z-coordinates of the gauss points
    AMP::LinearAlgebra::Vector::const_shared_ptr d_z_coord1;
    AMP::LinearAlgebra::Vector::const_shared_ptr d_z_coord2;

private:
    Discretization::createLibmeshElements libmeshElements;
};


} // Operator namespace
} // AMP namespace

#endif
