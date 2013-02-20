#include "ampmesh/Mesh.h"
#include "utils/Utilities.h"

#include "ampmesh/MultiMesh.h"
#include "ampmesh/SubsetMesh.h"
#include "ampmesh/structured/BoxMesh.h"
#ifdef USE_TRILINOS_STKMESH
    #include "ampmesh/STKmesh/STKMesh.h"
#endif
#ifdef USE_EXT_LIBMESH
    #include "ampmesh/libmesh/libMesh.h"
#endif
#ifdef USE_EXT_MOAB
    #include "ampmesh/moab/moabMesh.h"
#endif
#include "ampmesh/MeshElementVectorIterator.h"

#ifdef USE_AMP_VECTORS
    #include "vectors/Vector.h"
    #include "vectors/Variable.h"
    #include "vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
    #include "discretization/DOF_Manager.h"
    #include "discretization/simpleDOF_Manager.h"
#endif

#include <math.h>


namespace AMP {
namespace Mesh {

static unsigned int nextLocalMeshID = 1; 


/********************************************************
* Constructors                                          *
********************************************************/
Mesh::Mesh( const MeshParameters::shared_ptr &params_in )
{
    // Set the base properties
    AMP_ASSERT(sizeof(MeshElementID)==16);
    params = params_in;
    GeomDim = null;
    d_comm = params->comm;
    d_db = params->d_db;
    AMP_INSIST(d_comm!=AMP_MPI(AMP_COMM_NULL), "Communicator in mesh params must be non NULL ");
    setMeshID();
    d_name = "NULL";
}
Mesh::Mesh( const Mesh::shared_ptr &old_mesh )
{
    AMP_ERROR("Copy constructor is not Implimented Yet");
}


/********************************************************
* De-constructor                                        *
********************************************************/
Mesh::~Mesh()
{
}


/********************************************************
* Assignment operator                                   *
********************************************************/
Mesh Mesh::operator=(const Mesh& rhs)
{
    return rhs.copy();
}
Mesh Mesh::copy() const
{
    return Mesh(*this);
}


/********************************************************
* Create a mesh from the input database                 *
********************************************************/
boost::shared_ptr<AMP::Mesh::Mesh> Mesh::buildMesh( const MeshParameters::shared_ptr &params )
{
    boost::shared_ptr<AMP::Database> database = params->d_db;
    AMP_ASSERT(database!=NULL);
    AMP_INSIST(database->keyExists("MeshType"),"MeshType must exist in input database");
    AMP_INSIST(database->keyExists("MeshName"),"MeshName must exist in input database");
    std::string MeshType = database->getString("MeshType");
    std::string MeshName = database->getString("MeshName");
    boost::shared_ptr<AMP::Mesh::Mesh> mesh;
    if ( MeshType == std::string("Multimesh") ) {
        // The mesh is a multimesh
        mesh = boost::shared_ptr<AMP::Mesh::MultiMesh>(new AMP::Mesh::MultiMesh(params) );
    } else if ( MeshType == std::string("AMP") ) {
        // The mesh is a AMP mesh
        mesh = boost::shared_ptr<AMP::Mesh::BoxMesh>(new AMP::Mesh::BoxMesh(params) );
    } else if ( MeshType == std::string("libMesh") ) {
        // The mesh is a libmesh mesh
        #ifdef USE_EXT_LIBMESH
            mesh = boost::shared_ptr<AMP::Mesh::libMesh>(new AMP::Mesh::libMesh(params) );
        #else
            AMP_ERROR("AMP was compiled without support for libMesh");
        #endif
    } else if ( MeshType == std::string("STKMesh") ) {
        // The mesh is a libmesh mesh
        #ifdef USE_TRILINOS_STKMESH
            mesh = boost::shared_ptr<AMP::Mesh::STKMesh>(new AMP::Mesh::STKMesh(params) );
        #else
            AMP_ERROR("AMP was compiled without support for STKMesh");
        #endif
    } else if ( MeshType==std::string("moab") || MeshType==std::string("MOAB") ) {
        // The mesh is a MOAB mesh
        #ifdef USE_EXT_MOAB
            mesh = boost::shared_ptr<AMP::Mesh::moabMesh>(new AMP::Mesh::moabMesh(params) );
        #else
            AMP_ERROR("AMP was compiled without support for MOAB");
        #endif
    } else {
        // Unknown mesh type
        AMP_ERROR( std::string("Unknown mesh type (") + MeshType + std::string(")") );
    }
    mesh->setName(MeshName);
    return mesh;
}


/********************************************************
* Estimate the mesh size                                *
********************************************************/
size_t Mesh::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    boost::shared_ptr<AMP::Database> database = params->d_db;
    AMP_ASSERT(database!=NULL);
    size_t meshSize = 0;
    if ( database->keyExists("NumberOfElements") ) {
        // User specified the number of elements, this should override everything
        meshSize = (size_t) database->getInteger("NumberOfElements");
        // Adjust the number of elements by a weight if desired
        if ( database->keyExists("Weight") ) {
            double weight = database->getDouble("Weight");
            meshSize = (size_t) ceil(weight*((double)meshSize));
        }
        return meshSize;
    }
    // This is being called through the base class, call the appropriate function
    AMP_INSIST(database->keyExists("MeshType"),"MeshType must exist in input database");
    std::string MeshType = database->getString("MeshType");
    boost::shared_ptr<AMP::Mesh::Mesh> mesh;
    if ( MeshType == std::string("Multimesh") ) {
        // The mesh is a multimesh
        meshSize = AMP::Mesh::MultiMesh::estimateMeshSize(params);
    } else if ( MeshType == std::string("AMP") ) {
        // The mesh is a AMP mesh
        meshSize = AMP::Mesh::BoxMesh::estimateMeshSize(params);
    } else if ( MeshType == std::string("libMesh") ) {
        // The mesh is a libmesh mesh
        #ifdef USE_EXT_LIBMESH
            meshSize = AMP::Mesh::libMesh::estimateMeshSize(params);
        #else
            AMP_ERROR("AMP was compiled without support for libMesh");
        #endif
    } else if ( MeshType == std::string("STKMesh") ) {
        // The mesh is a stkMesh mesh
        #ifdef USE_TRILINOS_STKMESH
            meshSize = AMP::Mesh::STKMesh::estimateMeshSize(params);
        #else
            AMP_ERROR("AMP was compiled without support for STKMesh");
        #endif
    } else if ( database->keyExists("NumberOfElements") ) {
        int NumberOfElements = database->getInteger("NumberOfElements");
        meshSize = NumberOfElements;
    } else {
        // Unknown mesh type
        AMP_ERROR( "Unknown mesh type and NumberOfElements does not exist in database" );
    }
    return meshSize;
}


/********************************************************
* Function to set the mesh ID                           *
* This function will create a unique ID for every mesh. *
* To accomplish this goal, the ID will consist of the   *
* rank of the root processor (from the global comm),    *
* and the number of meshes created by that processor.   *
********************************************************/
void Mesh::setMeshID( )
{
    if ( d_comm.getRank()==0 ) {
        // Root will create the meshID
        AMP_MPI globalComm(AMP_COMM_WORLD);
        d_meshID = MeshID(globalComm.getRank(),nextLocalMeshID);
        nextLocalMeshID++;
    }
    // Broadcast the meshID to all processors
    d_meshID = d_comm.bcast(d_meshID,0);
}


/********************************************************
* Function to return the meshID composing the mesh      *
********************************************************/
std::vector<MeshID> Mesh::getAllMeshIDs() const
{
    return std::vector<MeshID>(1,d_meshID);
}
std::vector<MeshID> Mesh::getBaseMeshIDs() const
{
    return std::vector<MeshID>(1,d_meshID);
}
std::vector<MeshID> Mesh::getLocalMeshIDs() const
{
    return std::vector<MeshID>(1,d_meshID);
}
std::vector<MeshID> Mesh::getLocalBaseMeshIDs() const
{
    return std::vector<MeshID>(1,d_meshID);
}


/********************************************************
* Function to return the mesh with the given ID         *
********************************************************/
boost::shared_ptr<Mesh>  Mesh::Subset( MeshID meshID ) const {
    if ( d_meshID==meshID ) 
        return boost::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return boost::shared_ptr<Mesh>();
}


/********************************************************
* Function to return the mesh with the given name       *
********************************************************/
boost::shared_ptr<Mesh>  Mesh::Subset( std::string name ) const {
    if ( d_name==name ) 
        return boost::const_pointer_cast<Mesh>( shared_from_this() );
    else
        return boost::shared_ptr<Mesh>();
}


/********************************************************
* Function to subset a mesh using a mesh iterator       *
********************************************************/
boost::shared_ptr<Mesh> Mesh::Subset( const MeshIterator &iterator ) const
{
    boost::shared_ptr<const Mesh> this_mesh( shared_from_this() );
    boost::shared_ptr<SubsetMesh> mesh( new SubsetMesh( this_mesh, iterator ) );
    return mesh;
}


/********************************************************
* Function to return the element given an ID            *
********************************************************/
MeshElement Mesh::getElement ( const MeshElementID &elem_id ) const
{
    MeshID mesh_id = elem_id.meshID();
    AMP_INSIST(mesh_id==d_meshID,"mesh id must match the mesh id of the element");
    MeshIterator iterator = getIterator( elem_id.type() );
    for (size_t i=0; i<iterator.size(); i++) {
        if ( iterator->globalID() == elem_id )
            return *iterator;
        ++iterator;
    }
    return MeshElement();
}


/********************************************************
* Function to return parents of an element              *
********************************************************/
std::vector<MeshElement> Mesh::getElementParents ( const MeshElement&, const GeomType ) const
{
    AMP_ERROR("getElementParents is not implimented for the base class");
    return std::vector<MeshElement>();
}


/********************************************************
* Return the position vector                            *
********************************************************/
#ifdef USE_AMP_VECTORS
AMP::LinearAlgebra::Vector::shared_ptr  Mesh::getPositionVector( std::string name, const int gcw ) const
{
    #ifdef USE_AMP_DISCRETIZATION
        AMP::Discretization::DOFManager::shared_ptr DOFs = 
            AMP::Discretization::simpleDOFManager::create( 
            boost::const_pointer_cast<Mesh>(shared_from_this()), 
            AMP::Mesh::Vertex, gcw, PhysicalDim, true );
        AMP::LinearAlgebra::Variable::shared_ptr nodalVariable( new AMP::LinearAlgebra::Variable(name) );
        AMP::LinearAlgebra::Vector::shared_ptr position = AMP::LinearAlgebra::createVector( DOFs, nodalVariable, true );
        std::vector<size_t> dofs(PhysicalDim);
        AMP::Mesh::MeshIterator cur = DOFs->getIterator();
        AMP::Mesh::MeshIterator end = cur.end();
        while ( cur != end ) {
            AMP::Mesh::MeshElementID id = cur->globalID();
            std::vector<double> coord = cur->coord();
            DOFs->getDOFs( id, dofs );
            position->setValuesByGlobalID( dofs.size(), &dofs[0], &coord[0] );
            ++cur;
        }
        return position;
    #else
        AMP_ERROR("getPositionVector requires DISCRETIZATION");
        return AMP::LinearAlgebra::Vector::shared_ptr();
    #endif
}
#endif


/********************************************************
* Functions that aren't implimented for the base class  *
********************************************************/
boost::shared_ptr<Mesh> Mesh::Subset( Mesh & ) const
{
    AMP_ERROR("Subset is not implimented for the base class");
    return boost::shared_ptr<Mesh>();
}
MeshIterator Mesh::getIterator( const GeomType, const int ) const
{
    AMP_ERROR("getIterator is not implimented for the base class");
    return MeshIterator();
}
MeshIterator Mesh::getSurfaceIterator( const GeomType, const int ) const
{
    AMP_ERROR("getSurfaceIterator is not implimented for the base class");
    return MeshIterator();
}
std::vector<int> Mesh::getBoundaryIDs ( ) const
{
    AMP_ERROR("getBoundaryIDs is not implimented for the base class");
    return std::vector<int>();
}
MeshIterator Mesh::getBoundaryIDIterator ( const GeomType, const int, const int ) const
{
    AMP_ERROR("getBoundaryIDIterator is not implimented for the base class");
    return MeshIterator();
}
std::vector<int> Mesh::getBlockIDs ( ) const
{
    AMP_ERROR("getBlockIDs is not implimented for the base class");
    return std::vector<int>();
}
MeshIterator Mesh::getBlockIDIterator ( const GeomType, const int, const int ) const
{
    AMP_ERROR("getBlockIDIterator is not implimented for the base class");
    return MeshIterator();
}
size_t Mesh::numLocalElements( const GeomType type ) const
{
    AMP_ERROR("numLocalElements is not implimented for the base class");
    return 0;
}
size_t Mesh::numGlobalElements( const GeomType type ) const
{
    AMP_ERROR("numGlobalElements is not implimented for the base class");
    return 0;
}
size_t Mesh::numGhostElements( const GeomType type, int gcw ) const
{
    AMP_ERROR("numGhostElements is not implimented for the base class");
    return 0;
}
void Mesh::displaceMesh( std::vector<double> x )
{
    AMP_ERROR("displaceMesh is not implimented for the base class");
}
#ifdef USE_AMP_VECTORS
void Mesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr x )
{
    AMP_ERROR("displaceMesh is not implimented for the base class");
}
#endif


/********************************************************
* MeshIterator set operations                           *
********************************************************/
MeshIterator Mesh::getIterator( SetOP OP, const MeshIterator &A, const MeshIterator &B )
{
    if ( OP == Union ) {
        // Perform a union: A U B
        // Get the union using the mesh IDs
        std::set<MeshElementID> union_set;
        MeshIterator curA = A.begin();
        for (size_t i=0; i<A.size(); i++) {
            union_set.insert(curA->globalID());
            ++curA;
        }
        MeshIterator curB = B.begin();
        for (size_t i=0; i<B.size(); i++) {
            union_set.insert(curB->globalID());
            ++curB;
        }
        std::vector<MeshElementID> union_ids(union_set.begin(),union_set.end());
        // Create the iterator
        if ( union_ids.size()==A.size() ) {
            return MeshIterator(A);
        } else if ( union_ids.size()==B.size() ) {
            return MeshIterator(B);
        } else {
            boost::shared_ptr<std::vector<MeshElement> > elements( new std::vector<MeshElement>(union_ids.size()) );
            curA = A.begin();
            for (size_t i=0; i<A.size(); i++) {
                MeshElementID idA = curA->globalID();
                size_t index = Utilities::findfirst(union_ids,idA);
                if ( index == union_ids.size() ) { index--; }
                if ( union_ids[index] == idA )
                    (*elements)[index] = *curA;
                ++curA;
            }
            curB = B.begin();
            for (size_t i=0; i<B.size(); i++) {
                MeshElementID idB = curB->globalID();
                size_t index = Utilities::findfirst(union_ids,idB);
                if ( index == union_ids.size() ) { index--; }
                if ( union_ids[index] == idB )
                    (*elements)[index] = *curB;
                ++curB;
            }
            return MultiVectorIterator( elements, 0 );
        }
    } else if ( OP == Intersection ) {
        // Perform a intersection: A n B
        // Get the intersection using the mesh IDs
        if ( A.size()==0 || B.size()==0 )
            return MeshIterator();
        std::vector<MeshElementID> idA(A.size());
        MeshIterator curA = A.begin();
        for (size_t i=0; i<A.size(); i++) {
            idA[i] = curA->globalID();
            ++curA;
        }
        Utilities::quicksort(idA);
        std::vector<MeshElementID> intersection;
        intersection.reserve(B.size());
        MeshIterator curB = B.begin();
        for (size_t i=0; i<B.size(); i++) {
            MeshElementID idB = curB->globalID();
            size_t index = Utilities::findfirst(idA,idB);
            if ( index == idA.size() ) { index--; }
            if ( idA[index] == idB )
                intersection.push_back(idB);
            ++curB;
        }
        if ( intersection.empty() )
            return MeshIterator();
        // Sort the intersection and check for duplicates
        Utilities::quicksort(intersection);
        for (size_t i=1; i<intersection.size(); i++)
            AMP_ASSERT(intersection[i]!=intersection[i-1]);
        // Create the iterator
        if ( intersection.size()==A.size() ) {
            return MeshIterator(A);
        } else if ( intersection.size()==B.size() ) {
            return MeshIterator(B);
        } else {
            boost::shared_ptr<std::vector<MeshElement> > elements( new std::vector<MeshElement>(intersection.size()) );
            curB = B.begin();
            for (size_t i=0; i<B.size(); i++) {
                MeshElementID idB = curB->globalID();
                size_t index = Utilities::findfirst(intersection,idB);
                if ( index == intersection.size() ) { index--; }
                if ( intersection[index] == idB )
                    (*elements)[index] = *curB;
                ++curB;
            }
            return MultiVectorIterator( elements, 0 );
        }
    } else if ( OP == Complement ) {
        // Perform a Complement:  A - B
        // Get the compliment using the mesh IDs
        std::set<MeshElementID> compliment_set;
        MeshIterator curA = A.begin();
        for (size_t i=0; i<A.size(); i++) {
            compliment_set.insert(curA->globalID());
            ++curA;
        }
        MeshIterator curB = B.begin();
        for (size_t i=0; i<B.size(); i++) {
            compliment_set.erase(curB->globalID());
            ++curB;
        }
        std::vector<MeshElementID> compliment(compliment_set.begin(),compliment_set.end());
        if ( compliment.empty() )
            return MeshIterator();
        // Create the iterator
        if ( compliment.size()==A.size() ) {
            return MeshIterator(A);
        } else {
            boost::shared_ptr<std::vector<MeshElement> > elements( new std::vector<MeshElement>(compliment.size()) );
            curA = A.begin();
            for (size_t i=0; i<A.size(); i++) {
                MeshElementID idA = curA->globalID();
                size_t index = Utilities::findfirst(compliment,idA);
                if ( index == compliment.size() ) { index--; }
                if ( compliment[index] == idA )
                    (*elements)[index] = *curA;
                ++curA;
            }
            return MultiVectorIterator( elements, 0 );
        }
    } else {
        AMP_ERROR("Unknown set operation");
    }
    return MeshIterator();
}


} // Mesh namespace
} // AMP namespace

