#ifdef USE_AMP_MESH
#ifdef USE_AMP_VECTORS
#include "WriteSolutionToFile.h"

void printSolution( AMP::Mesh::Mesh::shared_ptr mesh,
                    AMP::LinearAlgebra::Vector::shared_ptr solVec,
                    std::string exeName )
{

    AMP::Discretization::DOFManager::shared_ptr dof_map = solVec->getDOFManager();

    AMP::Mesh::MeshIterator nd     = mesh->getIterator( AMP::Mesh::Vertex, 0 );
    AMP::Mesh::MeshIterator end_nd = nd.end();

    std::string fname = "results_" + exeName + ".txt";
    FILE *fp          = fopen( fname.c_str(), "w" );

    fprintf( fp, "%s\n\n", exeName.c_str() );
    fprintf( fp, "x, y, z,   u,  v,  w\n\n" );

    std::vector<size_t> dofs;
    for ( ; nd != end_nd; ++nd ) {
        std::vector<double> x = nd->coord();
        for ( size_t i = 0; i < x.size(); i++ ) fprintf( fp, "%lf, ", x[i] );
        fprintf( fp, ",    " );
        dof_map->getDOFs( nd->globalID(), dofs );
        for ( size_t i = 0; i < dofs.size(); i++ ) {
            double val = solVec->getLocalValueByGlobalID( dofs[i] );
            fprintf( fp, "%.13lf, ", val );
        } // end for i
        fprintf( fp, " \n" );
    } // end for nd

    fclose( fp );
}

#endif
#endif
