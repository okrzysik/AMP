#include "utils/AMP_MPI.h"
#include "operators/map/Map1Dto3D.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace Operator {


template <class T>
static T* getPtr( std::vector<T> &x ) {
    if ( x.size()== 0 )
        return NULL;
    return &x[0];
}


// Constructor
Map1Dto3D::Map1Dto3D(const boost::shared_ptr<OperatorParameters>& params):
    MapOperator (params)
{
    boost::shared_ptr<MapOperatorParameters> myparams = 
        boost::dynamic_pointer_cast<MapOperatorParameters>(params);
    d_MapMesh = myparams->d_MapMesh;
    reset(myparams);
}


void Map1Dto3D :: reset(const boost::shared_ptr<OperatorParameters>& params)
{
    boost::shared_ptr<MapOperatorParameters> myparams =
         boost::dynamic_pointer_cast<MapOperatorParameters>(params);

    AMP_INSIST( ((myparams.get()) != NULL), "NULL parameter" );
    AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );
    AMP_INSIST( !myparams->d_MapComm.isNull(), "NULL communicator" );
    d_MapComm = myparams->d_MapComm;
    d_MapMesh = myparams->d_MapMesh;

    computeZLocations();

    AMP_INSIST( myparams->d_db->keyExists("InputVariable"), "key not found" );
    std::string inpVar = myparams->d_db->getString("InputVariable");
    d_inpVariable.reset(new AMP::LinearAlgebra::Variable(inpVar));

    AMP_INSIST( myparams->d_db->keyExists("OutputVariable"), "key not found" );
    std::string outVar = myparams->d_db->getString("OutputVariable");
    d_outVariable.reset(new AMP::LinearAlgebra::Variable(outVar));

}


void Map1Dto3D::computeZLocations(){

    // Check that the mesh exists on some processors
    int N_mesh = d_MapComm.sumReduce<int>((d_MapMesh.get()!=NULL?1:0));
    AMP_ASSERT(N_mesh>0);
    
    // Get the local location of nodes on the boundary
    std::vector<double> t_zLocations;
    if ( d_MapMesh.get()!=NULL ) {
        // Get an iterator over the nodes on the boundary
        AMP::Mesh::MeshIterator  bnd = d_MapMesh->getIDsetIterator( AMP::Mesh::Vertex, d_boundaryId, 0 );
        AMP::Mesh::MeshIterator  end_bnd = bnd.end();

        double Xx=0;
        double Yy=0;
        if(bnd!=end_bnd) {
            std::vector<double> x = bnd->coord();
            AMP_ASSERT(x.size()==3);
            t_zLocations.push_back(x[2]);
            Xx = x[0];
            Yy = x[1];
            bnd++;
        }
        for( ; bnd != end_bnd; ++bnd) {
            std::vector<double> x = bnd->coord();
            if( (fabs(Xx-x[0]) <= 1.e-12) && (fabs(Yy-x[1]) <= 1.e-12) ){ 
                t_zLocations.push_back(x[2]);
            }
        }
    }

    // Make Z locations consistent across all processors.
    size_t myLen = t_zLocations.size();
    size_t totLen = d_MapComm.sumReduce(myLen);
    std::vector<double> zLocations( totLen );
    d_MapComm.allGather ( getPtr(t_zLocations), myLen , getPtr(zLocations) );

    // Add the coordinates (internally this will make sure the values are unique and sort)
    setZLocations( zLocations );

}


// Set the z locations
void Map1Dto3D::setZLocations( const std::vector<double> &z )
{
    const double TOL = 1e-12;
    d_zLocations = z;
    if ( d_zLocations.size() <= 1 )
        return;
    // Sort the entires
    AMP::Utilities::quicksort(d_zLocations);
    // Remove any duplicate entries
    size_t next=1;
    for (size_t i=1; i<d_zLocations.size(); i++) {
        if ( d_zLocations[i]-d_zLocations[next-1] > TOL ) {
            d_zLocations[next] = d_zLocations[i];
            next++;
        }
    }
    d_zLocations.resize(next);
}


void Map1Dto3D::apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &u,
     AMP::LinearAlgebra::Vector::shared_ptr  &r, const double , const double )
{ 

    if ( d_MapMesh.get()==NULL ) 
        return;

    AMP_ASSERT(u != NULL);
    AMP::LinearAlgebra::Vector::shared_ptr inputVec = u->subsetVectorForVariable(d_inpVariable);
    // AMP::LinearAlgebra::Vector::shared_ptr outputVec =  r->subsetVectorForVariable(d_outVariable);
    AMP_ASSERT(inputVec != NULL);
    AMP_ASSERT(outputVec != NULL);
    // outputVec->zero();

    std::vector<int> numFaceNodes(outputVec->getLocalSize(),0);

    const unsigned int numPoints = inputVec->getLocalSize();

    // Loop through the points on the surface
    AMP_ASSERT(d_zLocations.size()>=2);
    AMP_ASSERT(d_zLocations.size()==inputVec->getLocalSize());
    AMP_ASSERT(d_zLocations.size()==inputVec->getGlobalSize());
    const double TOL = 1e-12;
    AMP::Discretization::DOFManager::shared_ptr dof_map = outputVec->getDOFManager();
    std::vector<size_t> dofs(1);
    AMP::Mesh::MeshIterator  bnd = d_MapMesh->getIDsetIterator( AMP::Mesh::Vertex, d_boundaryId, 0 );
    const double z1 = d_zLocations[0]-TOL;
    const double z2 = d_zLocations[d_zLocations.size()-1]+TOL;
    for (size_t i=0; i<bnd.size(); i++) {
        dof_map->getDOFs( bnd->globalID(), dofs );
        std::vector<double> x = bnd->coord();
        AMP_INSIST(dofs.size()==1,"Map1Dto3D is currently implemented for scalar quantities only");
        AMP_INSIST(x.size()==3,"Map1Dto3D is currently implemented for 3D");
        // Perform linear interpolation
        double z = x[2];
        if ( z<z1 || z>z2 ) {
            // Point is outside interpolant, do nothing
            AMP_ERROR("Bad interpolant");
        }
        size_t index = AMP::Utilities::findfirst(d_zLocations,z);
        if ( index==0 ) { index = 1; }
        if ( index==d_zLocations.size() ) { index = d_zLocations.size()-1; }
        double dz = (z-d_zLocations[index-1])/(d_zLocations[index]-d_zLocations[index-1]);
        double f1 = inputVec->getValueByLocalID(index-1);
        double f2 = inputVec->getValueByLocalID(index);
        double dz2 = 1.0-dz;
        double f = dz2*f1 + dz*f2;
        outputVec->setValueByGlobalID(dofs[0],f);
        ++bnd;
    }

    if(d_iDebugPrintInfoLevel>4) {
        AMP::pout << "The input to Map1Dto3D " << std::endl;
        AMP::pout << inputVec << std::endl;
    }

    outputVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

    if(d_iDebugPrintInfoLevel>5) {
        AMP::pout << "The output to Map1Dto3D " << std::endl;
        AMP::pout << outputVec << std::endl;
    }

}


}
}

