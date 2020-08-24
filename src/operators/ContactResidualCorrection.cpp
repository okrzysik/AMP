
#include "ContactResidualCorrection.h"

namespace AMP {
namespace Operator {

void ContactResidualCorrection::apply( AMP::LinearAlgebra::Vector::const_shared_ptr,
                                       AMP::LinearAlgebra::Vector::shared_ptr r )
{
    AMP::LinearAlgebra::Vector::shared_ptr rMaster = r->subsetVectorForVariable( d_masterVariable );
    AMP::LinearAlgebra::Vector::shared_ptr rSlave  = r->subsetVectorForVariable( d_slaveVariable );

    AMP::Discretization::DOFManager::shared_ptr master_dof_map = rMaster->getDOFManager();
    AMP::Discretization::DOFManager::shared_ptr slave_dof_map  = rSlave->getDOFManager();

    for ( size_t i = 0; i < d_masterNodes.size(); i++ ) {
        std::vector<size_t> masterGlobalIds;
        std::vector<size_t> slaveGlobalIds;
        master_dof_map->getDOFs( d_masterNodes[i], masterGlobalIds );
        slave_dof_map->getDOFs( d_slaveNodes[i], slaveGlobalIds );
	const double zero = 0.0;
        for ( auto &elem : d_dofs[i] ) {
            double slaveVal = rSlave->getLocalValueByGlobalID( slaveGlobalIds[elem] );
            rMaster->addLocalValuesByGlobalID( 1, &masterGlobalIds[elem], &slaveVal );
            rSlave->setLocalValuesByGlobalID( 1, &slaveGlobalIds[elem], &zero );
            slaveVal = rSlave->getLocalValueByGlobalID( slaveGlobalIds[elem] );
        } // end for j
    }     // end for i
}
} // namespace Operator
} // namespace AMP
