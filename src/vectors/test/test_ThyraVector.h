#include "test_Vector.h"
#include "vectors/trilinos/ManagedThyraVector.h"
#include "vectors/trilinos/NativeThyraVector.h"
//#include "test_ThyraVectorTests.h"
#include "utils/AMP_MPI.h"


/// \cond UNDOCUMENTED

namespace AMP {
namespace unit_test {


class  NativeThyraFactory
{
public:
    typedef AMP::LinearAlgebra::ThyraVector                  vector;

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable() {
        return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::Variable ( "thyra" ) );
    }

    static AMP::LinearAlgebra::Vector::shared_ptr getVector() {
        AMP_MPI global_comm(AMP_COMM_WORLD);
        int local_size = 101;
        // Create an epetra vector
        #ifdef USE_EXT_MPI
            Epetra_MpiComm  comm = global_comm.getCommunicator();
        #else
            Epetra_SerialComm  comm;
        #endif
        RCP<Epetra_Map> epetra_map( new Epetra_Map(local_size,0,comm) );
        RCP<Epetra_Vector> epetra_v( new Epetra_Vector(epetra_map,true) );
        // Create a thyra vector from the epetra vector
        RCP<VectorSpaceBase<double> > space = create_VectorSpace( epetra_map );
        RCP<VectorBase<double> > thyra_v = create_Vector( epetra_v, space );
        // Create the NativeThyraVector
        NativeThyraVectorParameters params( thyra_v );
        boost::shared_ptr<NativeThyraVector> vec( new NativeThyraVector( params ) );
        return vec;
    }
};


}
}

/// \endcond
