#ifndef included_AMP_SubchannelToCladGPMap
#define included_AMP_SubchannelToCladGPMap

#include "discretization/createLibmeshElements.h"
#include "operators/map/SubchannelToCladMap.h"

namespace AMP {
namespace Operator {


typedef SubchannelToCladMapParameters SubchannelToCladGPMapParameters;


/**
 * \class  SubchannelToCladGPMap
 * \brief  A gauss-point version of SubchannelToCladMap
 */
class SubchannelToCladGPMap : public SubchannelToCladMap
{
public:
    /** \brief  Returns true if MapType = "SubchannelToCladGPMap"
     * \param[in] s  A string extracted from the MapType line in a MeshToMeshMap db
     * \return  True iff s == "SubchannelToCladMapParameters"
     */
    static bool validMapType( const std::string &s );

    /** \brief   Standard constructor
     * \param[in] params  Input parameters
     */
    explicit SubchannelToCladGPMap(
        const AMP::shared_ptr<AMP::Operator::OperatorParameters> &params );

    //! Destructor
    virtual ~SubchannelToCladGPMap();


protected:
    // For a given subchannel, fill the elements of interest using the coordinates
    virtual void fillReturnVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                                   double range[4],
                                   AMP::Mesh::Mesh::shared_ptr mesh,
                                   const std::vector<AMP::Mesh::MeshElementID> &ids,
                                   const std::vector<double> &z,
                                   const std::vector<double> &f );

private:
    Discretization::createLibmeshElements libmeshElements;

    struct gaussPointZCoord {
        double z[4];
    };
    std::vector<gaussPointZCoord>
    getGaussPoints( AMP::Mesh::Mesh::shared_ptr mesh,
                    const std::vector<AMP::Mesh::MeshElementID> &ids );
};


} // Operator namespace
} // AMP namespace

#endif
