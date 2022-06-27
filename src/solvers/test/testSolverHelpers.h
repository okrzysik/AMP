#include "AMP/vectors/Vector.h"

#include <fstream>
#include <functional>

// Check the solution of the form: T = a + b*z + c*z*z
inline bool checkAnalyticalSolution( const std::string &exeName,
                                     std::function<double( double, double, double )> fun,
                                     const AMP::Mesh::MeshIterator &iterator,
                                     std::shared_ptr<const AMP::LinearAlgebra::Vector> vec )
{
    // Serial execution
    bool passes = true;
    auto DOFmap = vec->getDOFManager();
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    for ( int i = 0; i < globalComm.getSize(); ++i ) {
        if ( globalComm.getRank() == i ) {
            std::string filename = "data_" + exeName;
            int rank             = globalComm.getRank();
            int nranks           = globalComm.getSize();
            auto omode           = std::ios_base::out;
            if ( rank > 0 )
                omode |= std::ios_base::app;
            std::ofstream file( filename.c_str(), omode );
            if ( rank == 0 ) {
                file << "(* x y z analytic calculated relative-error *)" << std::endl;
                file << "results={" << std::endl;
            }
            file.precision( 14 );

            size_t numNodes = 0, iNode = 0;
            for ( auto it = iterator.begin(); it != iterator.end(); ++it )
                numNodes++;

            double mse = 0.0;
            for ( auto it = iterator.begin(); it != iterator.end(); ++it ) {
                std::vector<size_t> gid;
                DOFmap->getDOFs( it->globalID(), gid );
                double cal = vec->getValueByGlobalID( gid[0] );
                double x   = ( it->coord() )[0];
                double y   = ( it->coord() )[1];
                double z   = ( it->coord() )[2];
                double sol = fun( x, y, z );
                double err =
                    fabs( cal - sol ) * 2. / ( cal + sol + std::numeric_limits<double>::epsilon() );
                mse += ( sol - cal ) * ( sol - cal );
                file << "{" << x << "," << y << "," << z << "," << sol << "," << cal << "," << err
                     << "}";
                if ( iNode < numNodes - 1 )
                    file << "," << std::endl;
                if ( fabs( cal - sol ) > cal * 1e-3 )
                    passes = false;
                iNode++;
            }

            if ( rank == nranks - 1 ) {
                file << "};" << std::endl;
                mse /= ( 1. * iNode );
                mse = sqrt( mse );
                file << "l2err = {" << iNode << "," << mse << "};\n";
            }
            file.close();
        }
        globalComm.barrier();
    }
    return passes;
}
