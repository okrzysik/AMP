#ifdef USE_AMP_MESH
    #ifdef USE_AMP_VECTORS
        #include "WriteSolutionToFile.h"

void printSolution( AMP::Mesh::Mesh::shared_ptr mesh,
                    AMP::LinearAlgebra::Vector::shared_ptr solVec,
                    std::string exeName )
{

    auto dof_map = solVec->getDOFManager();

    auto nd     = mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    auto end_nd = nd.end();

    std::string fname = "results_" + exeName + ".txt";
    FILE *fp          = fopen( fname.c_str(), "w" );

    fprintf( fp, "%s\n\n", exeName.c_str() );
    fprintf( fp, "x, y, z,   u,  v,  w\n\n" );

    std::vector<size_t> dofs;
    for ( ; nd != end_nd; ++nd ) {
        auto x = nd->coord();
        for ( auto &elem : x )
            fprintf( fp, "%lf, ", elem );
        fprintf( fp, ",    " );
        dof_map->getDOFs( nd->globalID(), dofs );
        for ( auto val : dofs ) {
            fprintf( fp, "%.13lf, ", static_cast<double>( val ) );
        } // end for i
        fprintf( fp, " \n" );
    } // end for nd

    fclose( fp );
}

    #endif
#endif
