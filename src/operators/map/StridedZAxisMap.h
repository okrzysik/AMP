#ifndef included_AMP_StridedZAxisMap
#define included_AMP_StridedZAxisMap

#include "operators/map/Map3to1to3.h"
#include "operators/map/Map3to1to3Parameters.h"
#include "operators/map/ScalarZAxisMap.h"


namespace AMP {
namespace Operator {

typedef AMP::Operator::Map3to1to3Parameters ScalarZAxisMapParameters;

/**
 * \class  StridedZAxisMap
 * \brief  A class used to reduce a 3D problem to 1D, transfer the solution, and map back to 3D
 * \details  This class inherites from ScalarZAxisMap, and performs a reduction from 3D to 1D,
 *    transfers the solution, then maps back to 3D.  It inherites entire
 *    functionality from ScalarZAxis except that the you need to do a subset
 *    and stride to pass the correct vector to the base class.
 */
class StridedZAxisMap : public AMP::Operator::ScalarZAxisMap
{
public:
    /** \brief   Standard constructor
     * \param[in] params  Input parameters
     */
    explicit StridedZAxisMap( const AMP::shared_ptr<AMP::Operator::OperatorParameters> &params );

    //! Destructor
    virtual ~StridedZAxisMap();

    /** \brief  Returns true if MapType = "ScalarZAxis"
     * \param[in] s  A string extracted from the MapType line in a MeshToMeshMap db
     * \return  True iff s == "ScalarZAxis"
     */
    static bool validMapType( const std::string &s );

    // Overload the apply
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f );

    virtual void setVector( AMP::LinearAlgebra::Vector::shared_ptr p );

protected:
    int d_inpDofs, d_inpStride, d_outDofs, d_outStride;
};


} // namespace Operator
} // namespace AMP

#endif
