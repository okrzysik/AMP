#include "ampmesh/libmesh/libMeshElement.h"
#include "utils/Utilities.h"

#include "boundary_info.h"

namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int libMeshElementTypeID = TYPE_HASH(libMeshElement);

// Functions to create new ids by mixing existing ids
static unsigned int mix2ids(unsigned int, unsigned int);
static unsigned int mix3ids(unsigned int, unsigned int, unsigned int);


/********************************************************
* Constructors                                          *
********************************************************/
libMeshElement::libMeshElement()
{
    typeID = libMeshElementTypeID;
    element = NULL;
    d_dim = -1;
    d_globalID = MeshElementID();
}
libMeshElement::libMeshElement(int dim, GeomType type, void* libmesh_element, 
    unsigned int rank, MeshID meshID, const libMesh* mesh)
{
    AMP_ASSERT(libmesh_element!=NULL);
    typeID = libMeshElementTypeID;
    element = NULL;
    d_dim = dim;
    d_rank = rank;
    d_mesh = mesh;
    d_meshID = meshID;
    ptr_element = libmesh_element;
    unsigned int local_id = (unsigned int)-1;
    unsigned int owner_rank = (unsigned int)-1;
    bool is_local=false;
    if ( type==Vertex ) {
        ::Node* node = (::Node*) ptr_element;
        local_id = node->id();
        owner_rank = node->processor_id();
        is_local = owner_rank==d_rank;
    } else if ( type==(GeomType) dim ) {
        ::Elem* elem = (::Elem*) ptr_element;
        AMP_ASSERT(elem->n_neighbors()<100);
        local_id = elem->id();
        owner_rank = elem->processor_id();
        is_local = owner_rank==d_rank;
    } else {
        AMP_ERROR("Unreconized element");
    }
    d_globalID = MeshElementID(is_local,type,local_id,owner_rank,meshID);
}
libMeshElement::libMeshElement(int dim, GeomType type, boost::shared_ptr< ::Elem > libmesh_element, 
    unsigned int rank, MeshID meshID, const libMesh* mesh)
{
    AMP_ASSERT(libmesh_element.get()!=NULL);
    typeID = libMeshElementTypeID;
    element = NULL;
    d_dim = dim;
    d_rank = rank;
    d_mesh = mesh;
    d_meshID = meshID;
    ptr2 = libmesh_element;
    ptr_element = libmesh_element.get();
    unsigned int local_id = (unsigned int)-1;
    unsigned int owner_rank = (unsigned int)-1;
    bool is_local=false;
    if ( type==Vertex ) {
        ::Node* node = (::Node*) ptr_element;
        local_id = node->id();
        owner_rank = node->processor_id();
        is_local = owner_rank==d_rank;
    } else {
        ::Elem* elem = (::Elem*) ptr_element;
        local_id = elem->id();
        owner_rank = elem->processor_id();
        is_local = owner_rank==d_rank;
    }
    d_globalID = MeshElementID(is_local,type,local_id,owner_rank,meshID);
}
libMeshElement::libMeshElement(const libMeshElement& rhs)
{
    typeID = libMeshElementTypeID;
    element = NULL;
    d_globalID = rhs.d_globalID;
    d_dim = rhs.d_dim;
    ptr_element = rhs.ptr_element;
    ptr2 = rhs.ptr2;
    d_rank = rhs.d_rank;
    d_mesh = rhs.d_mesh;
    d_meshID = rhs.d_meshID;
}
libMeshElement& libMeshElement::operator=(const libMeshElement& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    this->typeID = libMeshElementTypeID;
    this->element = NULL;
    this->d_globalID = rhs.d_globalID;
    this->d_dim = rhs.d_dim;
    this->ptr_element = rhs.ptr_element;
    this->ptr2 = rhs.ptr2;
    this->d_rank = rhs.d_rank;
    this->d_mesh = rhs.d_mesh;
    this->d_meshID = rhs.d_meshID;
    return *this;
}


/****************************************************************
* De-constructor                                                *
****************************************************************/
libMeshElement::~libMeshElement()
{
    element = NULL;
}


/****************************************************************
* Function to clone the element                                 *
****************************************************************/
MeshElement* libMeshElement::clone() const
{
    return new libMeshElement(*this);
}


/****************************************************************
* Function to get the elements composing the current element    *
****************************************************************/
std::vector<MeshElement> libMeshElement::getElements(const GeomType type) const
{
    AMP_INSIST(type<=d_globalID.type(),"sub-elements must be of a smaller or equivalent type");
    std::vector<MeshElement> children(0);
    ::Elem* elem = (::Elem*) ptr_element;
    if ( d_globalID.type()==Vertex ) {
        // A vertex does not have children, return itself
        if ( type!=Vertex )
            AMP_ERROR("A vertex is the base element and cannot have and sub-elements");
        children.resize(1);
        children[0] = *this;
    } else if ( type==d_globalID.type() ) {
        // Return the children of the current element
        if ( elem->has_children() ) {
            children.resize(elem->n_children());
            for (unsigned int i=0; i<children.size(); i++)
                children[i] = libMeshElement( d_dim, type, (void*)elem->child(i), d_rank, d_meshID, d_mesh );
        } else {
            children.resize(1);
            children[0] = *this;
        }
    } else if ( type==Vertex ) {
        // Return the nodes of the current element
        children.resize(elem->n_nodes());
        for (unsigned int i=0; i<children.size(); i++)
            children[i] = libMeshElement( d_dim, type, (void*)elem->get_node(i), d_rank, d_meshID, d_mesh );
    } else {
        // Return the children
        if ( type==Edge )
            children.resize(elem->n_edges());
        else if ( type==Face )
            children.resize(elem->n_faces());
        else
            AMP_ERROR("Internal error");
        for (unsigned int i=0; i<children.size(); i++) {
            // We need to build a valid element
            ::AutoPtr< ::Elem > tmp;
            if ( type==Edge )
                tmp = elem->build_edge(i);
            else if ( type==Face )
                tmp = elem->build_side(i,false);
            else
                AMP_ERROR("Internal error");
            boost::shared_ptr< ::Elem > element( tmp.release() );
            // We need to generate a vaild id and owning processor
            unsigned int n_node_min = ((unsigned int)type) + 1;
            AMP_ASSERT(element->n_nodes()>=n_node_min);
            std::vector< ::Node* > nodes(element->n_nodes());
            std::vector<unsigned int> node_ids(element->n_nodes());
            for (size_t j=0; j<nodes.size(); j++) {
                nodes[j] = element->get_node(j);
                node_ids[j] = nodes[j]->id();
            }
            AMP::Utilities::quicksort(node_ids,nodes);
            element->processor_id() = nodes[0]->processor_id();
            if ( n_node_min==2 ) 
                element->set_id() = mix2ids(node_ids[0],node_ids[1]);
            else if ( n_node_min==3 ) 
                element->set_id() = mix3ids(node_ids[0],node_ids[1],node_ids[2]);
            else
                AMP_ERROR("Internal error");
            // Create the libMeshElement
            children[i] = libMeshElement( d_dim, type, element, d_rank, d_meshID, d_mesh );            
        }
    }
    return children;
}


/****************************************************************
* Function to get the neighboring elements                      *
****************************************************************/
std::vector< MeshElement::shared_ptr > libMeshElement::getNeighbors() const
{
    std::vector< MeshElement::shared_ptr > neighbors(0);
    if ( d_globalID.type()==Vertex ) {
        // Return the neighbors of the current node
        std::vector< ::Node* > neighbor_nodes = d_mesh->getNeighborNodes( d_globalID );
        neighbors.resize(neighbor_nodes.size(),MeshElement::shared_ptr());
        for (size_t i=0; i<neighbor_nodes.size(); i++) {
            // There are no NULL neighbors
            boost::shared_ptr<libMeshElement> neighbor(new libMeshElement( d_dim, Vertex, (void*)neighbor_nodes[i], d_rank, d_meshID, d_mesh ) );
            neighbors[i] = neighbor;
        }
    } else if ( (int) d_globalID.type() == d_dim ) {
        // Return the neighbors of the current element
        ::Elem* elem = (::Elem*) ptr_element;
        //if ( elem->n_neighbors()==0 )
        //    AMP_ERROR("Element has not neighbors, this could indicate a problem with the mesh");
    	neighbors.resize(elem->n_neighbors());
        for (size_t i=0; i<neighbors.size(); i++) {
            void *neighbor_elem = (void*) elem->neighbor(i);
            boost::shared_ptr<libMeshElement> neighbor;
            if ( neighbor_elem != NULL )
                neighbor = boost::shared_ptr<libMeshElement>(new libMeshElement( d_dim, d_globalID.type(), neighbor_elem, d_rank, d_meshID, d_mesh ) );
            neighbors[i] = neighbor;
        }
    } else {
        // We constructed a temporary libmesh object and do not have access to the neighbor info
    }
    return neighbors;
}


/****************************************************************
* Functions to get basic element properties                     *
****************************************************************/
double libMeshElement::volume() const
{
    if ( d_globalID.type() == Vertex )
        AMP_ERROR("volume is is not defined Nodes");
    ::Elem* elem = (::Elem*) ptr_element;
    return elem->volume();
}
std::vector<double> libMeshElement::coord() const
{
    if ( d_globalID.type() != Vertex )
        AMP_ERROR("coord is only defined for Nodes");
    ::Node* node = (::Node*) ptr_element;
    std::vector<double> x(d_dim,0.0);
    for (int i=0; i<d_dim; i++)
        x[i] = (*node)(i);
    return x;
}
std::vector<double> libMeshElement::centroid() const
{
    if ( d_globalID.type()==Vertex )
        return coord();
    ::Elem* elem = (::Elem*) ptr_element;
    ::Point center = elem->centroid();
    std::vector<double> x(d_dim,0.0);
    for (int i=0; i<d_dim; i++)
        x[i] = center(i);
    return x;
}
bool libMeshElement::containsPoint( const std::vector<double> &pos, double TOL ) const
{
    if ( d_globalID.type()==Vertex ) {
        //double dist = 0.0;
        std::vector<double> point = this->coord();
        double dist2 = 0.0;
        for (size_t i=0; i<point.size(); i++)
            dist2 += (point[i]-pos[i])*(point[i]-pos[i]);
        return dist2<=TOL*TOL;
    }
    ::Elem* elem = (::Elem*) ptr_element;
    ::Point point(pos[0],pos[1],pos[2]);
    return elem->contains_point(point,TOL);
}
bool libMeshElement::isOnSurface() const
{
    GeomType type = d_globalID.type();
    MeshElement search = MeshElement(*this);
    if ( d_globalID.is_local() ) {
        const std::vector<MeshElement> &data = *(d_mesh->d_localSurfaceElements[type]);
        if ( data.empty() )
            return false;   // There are no elements on the surface for this processor
        size_t index = AMP::Utilities::findfirst( data, search );
        if ( index < data.size() ) {
            if ( d_mesh->d_localSurfaceElements[type]->operator[](index).globalID() == d_globalID )
                return true;
        }
    } else {
        const std::vector<MeshElement> &data = *(d_mesh->d_ghostSurfaceElements[type]);
        if ( data.empty() )
            return false;   // There are no elements on the surface for this processor
        size_t index = AMP::Utilities::findfirst( data, search );
        if ( index < data.size() ) {
            if ( d_mesh->d_ghostSurfaceElements[type]->operator[](index).globalID() == d_globalID )
                return true;
        }
    }
    return false;
}
bool libMeshElement::isOnBoundary(int id) const
{
    GeomType type = d_globalID.type();
    bool on_boundary = false;
    boost::shared_ptr< ::Mesh> d_libMesh = d_mesh->getlibMesh();
    if ( type==Vertex ) {
        // Entity is a libmesh node
        ::Node* node = (::Node*) ptr_element;
        std::vector< short int > bids = d_libMesh->boundary_info->boundary_ids(node);
        for (size_t i=0; i<bids.size(); i++) {
            if ( bids[i]==id )
                on_boundary = true;
        }
    } else if ( (int)type==d_dim ) {
        // Entity is a libmesh node
        ::Elem* elem = (::Elem*) ptr_element;
        unsigned int side = d_libMesh->boundary_info->side_with_boundary_id(elem,id);
        if ( side != static_cast<unsigned int>(-1) )
            on_boundary = true;
    } else  {
        // All other entities are on the boundary iff all of their verticies are on the surface
        std::vector<MeshElement> nodes = this->getElements(Vertex);
        on_boundary = true;
        for (size_t i=0; i<nodes.size(); i++) 
            on_boundary = on_boundary && nodes[i].isOnBoundary(id);
    }
    return on_boundary;
}
bool libMeshElement::isInBlock(int id) const
{
   GeomType type = d_globalID.type();
   bool in_block = false;
    if ( type==Vertex ) {
        // Entity is a libmesh node
        AMP_ERROR("isInBlock is not currently implimented for anything but elements");
    } else if ( (int)type==d_dim ) {
        // Entity is a libmesh node
        ::Elem* elem = (::Elem*) ptr_element;
        in_block = elem->subdomain_id() == id;
    } else  {
        // All other entities are on the boundary iff all of their verticies are on the surface
        AMP_ERROR("isInBlock is not currently implimented for anything but elements");
    }
    return in_block;
}


/****************************************************************
* Functions to mix ids                                          *
****************************************************************/
unsigned int mix2ids(unsigned int a, unsigned int b)
{
    return mix3ids(a,b,0);
    //return (a<<15) + b;
}
unsigned int mix3ids(unsigned int a, unsigned int b, unsigned int c)
{
    //return (a<<22) + (b<<11)  + c;
    a=a-b;  a=a-c;  a=a^(c>>13);
    b=b-c;  b=b-a;  b=b^(a<<8);
    c=c-a;  c=c-b;  c=c^(b>>13);
    a=a-b;  a=a-c;  a=a^(c>>12);
    b=b-c;  b=b-a;  b=b^(a<<16);
    c=c-a;  c=c-b;  c=c^(b>>5);
    a=a-b;  a=a-c;  a=a^(c>>3);
    b=b-c;  b=b-a;  b=b^(a<<10);
    c=c-a;  c=c-b;  c=c^(b>>15);
    return c;
}


} // Mesh namespace
} // AMP namespace

