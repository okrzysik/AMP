#include "ampmesh/MultiMesh.h"
#include "ampmesh/MultiIterator.h"
#include "ampmesh/MeshElement.h"
#include "utils/Utilities.h"
#include "utils/Database.h"
#include "utils/MemoryDatabase.h"


namespace AMP {
namespace Mesh {


// Function to check if a string contains another string as the prefix
static bool check_prefix(std::string prefix, std::string str) {
    if ( str.size() < prefix.size() )
        return false;
    if ( str.compare(0,prefix.size(),prefix)==0 )
        return true;
    return false;
}


// Some template functions for vector math
template <class T> 
static inline T vecmin(std::vector<T> x) {
    T min = x[0];
    for (size_t i=1; i<x.size(); i++)
        min = (x[i]<min) ? x[i] : min;
    return min;
}
template <class T> 
static inline T vecmax(std::vector<T> x) {
    T max = x[0];
    for (size_t i=1; i<x.size(); i++)
        max = (x[i]>max) ? x[i] : max;
    return max;
}


// Misc function declerations
static void copyKey(boost::shared_ptr<AMP::Database>&,boost::shared_ptr<AMP::Database>&,std::string,int,int);


/********************************************************
* Constructors                                          *
********************************************************/
MultiMesh::MultiMesh( const MeshParameters::shared_ptr &params_in ):
    Mesh(params_in)
{
    // Create a database for each mesh within the multimesh
    AMP_ASSERT(d_db!=NULL);
    std::vector<boost::shared_ptr<AMP::Database> > meshDatabases = MultiMesh::createDatabases(d_db);
    // Create the input parameters and get the approximate number of elements for each mesh
    std::vector<size_t> meshSizes(meshDatabases.size(),0);
    std::vector< boost::shared_ptr<AMP::Mesh::MeshParameters> > meshParameters(meshDatabases.size());
    for (size_t i=0; i<meshDatabases.size(); i++) {
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(meshDatabases[i]));
        params->setComm(AMP::AMP_MPI(AMP_COMM_SELF));
        meshParameters[i] = params;
        meshSizes[i] = AMP::Mesh::Mesh::estimateMeshSize(params);
    }
    size_t totalMeshSize = 0;
    for (size_t i=0; i<meshSizes.size(); i++)
        totalMeshSize += meshSizes[i];
    // Determine the load balancing we want to use and create the communicator for the parameters for each mesh
    int method = 1;
    if ( method == 1 ) {
        // This is the simplist load balancing, all meshes exist on the same communicator
        // This is equivalant to the old DomainDecomposition=0
        for (size_t i=0; i<meshParameters.size(); i++)
            meshParameters[i]->setComm(comm);
    } else if ( method == 2 ) {
        // This tries not mak each comm as small as possible while balancing the elements.  
        // No communicators may overlap, and only a simple load balancing scheme is used.
        // This is similar to the old DomainDecomposition=1, but does NOT guarantee that a
        // processor contains only one mesh.
        // First get the ~ number of processors needed for each mesh
        std::vector<double> N_procs(meshParameters.size(),0);
        for (size_t i=0; i<meshSizes.size(); i++)
            N_procs[i] = ((double)meshSizes[i])/((double)totalMeshSize)*((double)comm.getSize());
        // Group the meshes into communicator groups
        std::vector<int> group(meshParameters.size(),-1);
        std::vector<int> commSize(0);
        if ( vecmin(N_procs)>=0.5 && comm.getSize()>=(int)meshSizes.size() ) {
            // We have enough processors for each mesh to have a seperate group 
            AMP_ERROR("Not finished");
        } else {
            // We need to combine some processors into larger groups for their communicators
            AMP_ERROR("Not finished");
        }
        AMP_ERROR("Not finished");
    } else {
        AMP_ERROR("Unknown load balancer");
    }
    // Create the meshes
    d_meshes = std::vector<AMP::Mesh::Mesh::shared_ptr>(meshParameters.size());
    for (size_t i=0; i<meshParameters.size(); i++)
        d_meshes[i] = AMP::Mesh::Mesh::buildMesh(meshParameters[i]);
}


/********************************************************
* De-constructor                                        *
********************************************************/
MultiMesh::~MultiMesh()
{
}


/********************************************************
* Function to copy the mesh                             *
********************************************************/
Mesh MultiMesh::copy() const
{
    return MultiMesh(*this);
}


/********************************************************
* Function to estimate the mesh size                    *
********************************************************/
size_t MultiMesh::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    // Create a database for each mesh within the multimesh
    boost::shared_ptr<AMP::Database> database = params->getDatabase();
    AMP_ASSERT(database.get()!=NULL);
    std::vector<boost::shared_ptr<AMP::Database> > meshDatabases = MultiMesh::createDatabases(database);
    // Create the input parameters and get the approximate number of elements for each mesh
    size_t totalMeshSize = 0;
    std::vector< boost::shared_ptr<AMP::Mesh::MeshParameters> > meshParameters(meshDatabases.size());
    for (size_t i=0; i<meshDatabases.size(); i++) {
        boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(meshDatabases[i]));
        params->setComm(AMP::AMP_MPI(AMP_COMM_SELF));
        meshParameters.push_back(params);
        totalMeshSize += AMP::Mesh::Mesh::estimateMeshSize(params);
    }
    return totalMeshSize;
}


/********************************************************
* Function to create the databases for the meshes       *
* within the multimesh.                                 *
********************************************************/
std::vector<boost::shared_ptr<AMP::Database> >  MultiMesh::createDatabases(boost::shared_ptr<AMP::Database> database)
{
    // Find all of the meshes in the database
    AMP_ASSERT(database!=NULL);
    AMP_INSIST(database->keyExists("MeshDatabasePrefix"),"MeshDatabasePrefix must exist in input database");
    AMP_INSIST(database->keyExists("MeshArrayDatabasePrefix"),"MeshArrayDatabasePrefix must exist in input database");
    std::string MeshPrefix = database->getString("MeshDatabasePrefix");
    std::string MeshArrayPrefix = database->getString("MeshArrayDatabasePrefix");
    AMP_ASSERT(!check_prefix(MeshPrefix,MeshArrayPrefix));
    AMP_ASSERT(!check_prefix(MeshArrayPrefix,MeshPrefix));
    std::vector<std::string> keys = database->getAllKeys();
    std::vector<std::string> meshes, meshArrays;
    for (size_t i=0; i<keys.size(); i++) {
        if ( check_prefix(MeshPrefix,keys[i]) ) {
            meshes.push_back(keys[i]);
        } else if ( check_prefix(MeshArrayPrefix,keys[i]) ) {
            meshArrays.push_back(keys[i]);
        }
    }
    // Create the basic databases for each mesh
    std::vector<boost::shared_ptr<AMP::Database> > meshDatabases;
    for (size_t i=0; i<meshes.size(); i++) {
        // We are dealing with a single mesh object (it might be a multimesh), use the existing database
        boost::shared_ptr<AMP::Database> database2 = database->getDatabase(meshes[i]);
        meshDatabases.push_back(database2);
    }
    for (size_t i=0; i<meshArrays.size(); i++) {
        // We are dealing with an array of meshes, create a database for each
        boost::shared_ptr<AMP::Database> database1 = database->getDatabase(meshArrays[i]);
        int N = database1->getInteger("N");
        for (int j=0; j<N; j++) {
            // Create a new database from the existing database
            boost::shared_ptr<AMP::Database> database2( new AMP::MemoryDatabase(meshArrays[i]) );
            keys = database1->getAllKeys();
            for (size_t k=0; k<keys.size(); k++) {
                if ( keys[k].compare("N")==0 || keys[k].compare("iterator")==0 || keys[k].compare("indicies")==0 ) {
                    // These keys are used by the mesh-array and should not be copied
                } else {
                    // We need to copy the key
                    copyKey(database1,database2,keys[k],N,i);
                }
            }
            // We still need to perform the text search and replace

            boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database2));
            meshDatabases.push_back(database2);
        }
    }
    return meshDatabases;
}


/********************************************************
* Return basic mesh info                                *
********************************************************/
size_t MultiMesh::numLocalElements( const GeomType type ) const
{
    // Should we cache this?
    size_t N = 0;
    for (size_t i=0; i<d_meshes.size(); i++)
        N += d_meshes[i]->numLocalElements(type);
    return N;
}
size_t MultiMesh::numGlobalElements( const GeomType type ) const
{
    // Should we cache this?
    size_t N = 0;
    for (size_t i=0; i<d_meshes.size(); i++)
        N += d_meshes[i]->numGlobalElements(type);
    return N;
}
size_t MultiMesh::numGhostElements( const GeomType type, int gcw ) const
{
    // Should we cache this?
    size_t N = 0;
    for (size_t i=0; i<d_meshes.size(); i++)
        N += d_meshes[i]->numGhostElements(type,gcw);
    return N;
}
GeomType MultiMesh::getGeomType( ) const
{
    // Should we cache this?
    GeomType type_max = null;
    for (size_t i=0; i<d_meshes.size(); i++) {
        GeomType type = d_meshes[i]->getGeomType();
        type_max = (type>type_max) ? type : type_max;
    }
    return type_max;
}


/********************************************************
* Return mesh iterators                                 *
********************************************************/
MeshIterator MultiMesh::getIterator( const GeomType type, const int gcw )
{
    std::vector<boost::shared_ptr<MeshIterator> > iterators(d_meshes.size());
    for (size_t i=0; i<d_meshes.size(); i++) {
        boost::shared_ptr<MeshIterator> iterator_ptr( new MeshIterator(d_meshes[i]->getIterator(type,gcw)) );
        iterators[i] = iterator_ptr;
    }
    return MultiIterator(iterators);
}
std::vector<int> MultiMesh::getIDSets ( )
{
    AMP_ERROR("Not Implimented Yet");
    return std::vector<int>();
}
MeshIterator MultiMesh::getIDsetIterator ( const GeomType, const int, const int )
{
    AMP_ERROR("Not Implimented Yet");
    return MeshIterator();
}


/********************************************************
* Function to copy a key from database 1 to database 2  *
* If the key is an array of size N, it will only copy   *
* the ith value.                                        *
********************************************************/
static void copyKey(boost::shared_ptr<AMP::Database> &database1, 
    boost::shared_ptr<AMP::Database> &database2, std::string key, int N, int i ) 
{
    int size = database1->getArraySize(key);
    AMP::Database::DataType type = database1->getArrayType(key);
    switch (type) {
        case AMP::Database::AMP_INVALID: {
            // We don't know what this is
            AMP_ERROR("Invalid database object");
            } break;
        case AMP::Database::AMP_DATABASE: {
            // Copy the database
            AMP_ERROR("Not programmed for databases yet");
            } break;
        case AMP::Database::AMP_BOOL: {
            // Copy a bool
            std::vector<unsigned char> data = database1->getBoolArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putBool(key,data[i]);
            else
                database2->putBoolArray(key,data);
            } break;
        case AMP::Database::AMP_CHAR: {
            // Copy a char
            std::vector<char> data = database1->getCharArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size ) {
                database2->putChar(key,data[i]);
            } else {
                // We need to try a search and replace
                database2->putCharArray(key,data);
            }
            } break;
        case AMP::Database::AMP_INT: {
            // Copy an int
            std::vector<int> data = database1->getIntegerArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putInteger(key,data[i]);
            else
                database2->putIntegerArray(key,data);
            } break;
        case AMP::Database::AMP_COMPLEX: {
            // Copy a complex number
            std::vector< std::complex<double> > data = database1->getComplexArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putComplex(key,data[i]);
            else
                database2->putComplexArray(key,data);
            } break;
        case AMP::Database::AMP_DOUBLE: {
            // Copy a double number
            std::vector<double> data = database1->getDoubleArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putDouble(key,data[i]);
            else
                database2->putDoubleArray(key,data);
            } break;
        case AMP::Database::AMP_FLOAT: {
            // Copy a float
            std::vector<float> data = database1->getFloatArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putFloat(key,data[i]);
            else
                database2->putFloatArray(key,data);
            } break;
        case AMP::Database::AMP_STRING: {
            // Copy a string
            std::vector<std::string> data = database1->getStringArray(key);
            AMP_INSIST((int)data.size()==size,"Array size does not match key size");
            if ( N == size )
                database2->putString(key,data[i]);
            else
                database2->putStringArray(key,data);
            } break;
        case AMP::Database::AMP_BOX: {
            // Copy a box
            AMP_ERROR("Not programmed for boxes yet");
            } break;
        default:
            AMP_ERROR("Unknown key type");
    }
}


} // Mesh namespace
} // AMP namespace

