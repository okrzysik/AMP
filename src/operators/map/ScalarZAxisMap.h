#ifndef included_AMP_ScalarZAxisMap
#define included_AMP_ScalarZAxisMap

#include "AMP/operators/map/Map3to1to3.h"
#include "AMP/operators/map/Map3to1to3Parameters.h"


namespace AMP {
namespace Operator {


typedef AMP::Operator::Map3to1to3Parameters ScalarZAxisMapParameters;


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
    static bool validMapType( const std::string &s );

    /** \brief  Typedef to identify the parameters class of this operator
     */
    typedef ScalarZAxisMapParameters Parameters;

    //!  The base tag used in communication.
    static const int CommTagBase = 20000;

    /** \brief   Standard constructor
     * \param[in] params  Input parameters
     */
    explicit ScalarZAxisMap( const std::shared_ptr<AMP::Operator::OperatorParameters> &params );

    //! Destructor
    virtual ~ScalarZAxisMap();

    //! Return the name of the operator
    std::string type() const override { return "ScalarZAxisMap"; }

protected:
    // Implemented buildMap routine
    virtual std::multimap<double, double> buildMap( AMP::LinearAlgebra::Vector::const_shared_ptr,
                                                    const AMP::Mesh::Mesh::shared_ptr,
                                                    const AMP::Mesh::MeshIterator & ) override;

    // Implimented buildReturn routine
    virtual void buildReturn( AMP::LinearAlgebra::Vector::shared_ptr,
                              const AMP::Mesh::Mesh::shared_ptr,
                              const AMP::Mesh::MeshIterator &,
                              const std::map<double, double> & ) override;
};


} // namespace Operator
} // namespace AMP

#endif
