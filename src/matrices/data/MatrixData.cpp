#include "AMP/matrices/data/MatrixData.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/AMPManager.h"

namespace AMP::LinearAlgebra {


/********************************************************
 * Constructors/Destructor                               *
 ********************************************************/
MatrixData::MatrixData() { AMPManager::incrementResource( "MatrixData" ); }
MatrixData::MatrixData( std::shared_ptr<MatrixParametersBase> params ) : d_pParameters( params )
{
    AMPManager::incrementResource( "MatrixData" );
}
MatrixData::~MatrixData() { AMPManager::decrementResource( "MatrixData" ); }


/********************************************************
 * Get the number of rows/columns in the matrix          *
 ********************************************************/
size_t MatrixData::numLocalRows() const
{
    auto DOF = getLeftDOFManager();
    return DOF->numLocalDOF();
}
size_t MatrixData::numGlobalRows() const
{
    auto DOF = getLeftDOFManager();
    return DOF->numGlobalDOF();
}
size_t MatrixData::numLocalColumns() const
{
    auto DOF = getRightDOFManager();
    return DOF->numLocalDOF();
}
size_t MatrixData::numGlobalColumns() const
{
    auto DOF = getRightDOFManager();
    return DOF->numGlobalDOF();
}


/********************************************************
 * Get iterators                                         *
 ********************************************************/
size_t MatrixData::beginRow() const
{
    auto DOF = getLeftDOFManager();
    return DOF->beginDOF();
}

size_t MatrixData::endRow() const
{
    auto DOF = getLeftDOFManager();
    return DOF->endDOF();
}

size_t MatrixData::beginCol() const
{
    auto DOF = getRightDOFManager();
    return DOF->beginDOF();
}
size_t MatrixData::endCol() const
{
    auto DOF = getRightDOFManager();
    return DOF->endDOF();
}

} // namespace AMP::LinearAlgebra
