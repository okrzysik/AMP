#ifndef included_AMP_NodeToNodeMapParameters
#define included_AMP_NodeToNodeMapParameters

#include "ampmesh/Mesh.h"
#include "operators/map/AsyncMapOperatorParameters.h"

namespace AMP {
namespace Operator {

class NodeToNodeMapParameters : public AMP::Operator::AsyncMapOperatorParameters
{
public:
    std::vector<AMP::Mesh::MeshElementID> d_ids;
    std::vector<double> d_disps;
    std::map<size_t, size_t> d_RemoteToLocalId;
    size_t d_NumPartners;


    explicit NodeToNodeMapParameters( const AMP::shared_ptr<AMP::Database> &db )
        : AsyncMapOperatorParameters( db ), d_NumPartners( 0 )
    {
    }
};
} // namespace Operator
} // namespace AMP

#endif
