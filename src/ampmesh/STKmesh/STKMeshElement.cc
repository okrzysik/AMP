#include "ampmesh/STKmesh/STKMeshElement.h"
#include "utils/Utilities.h"

#include "stk_mesh/fem/FEMHelpers.hpp"
#include "stk_mesh/base/FieldData.hpp"
#include "stk_mesh/fem/CellTopology.hpp"
#include "Teuchos_RCP.hpp"

#include "Intrepid_Types.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"


namespace AMP {
namespace Mesh {


// Create a unique id for this class
namespace {
unsigned int STKMeshElementTypeID() {
    static const unsigned int id = TYPE_HASH(STKMeshElement);
    return id;
}
AMP::Mesh::GeomType geom_type(const stk::mesh::EntityRank rank) {
    switch (rank) {
      case stk::mesh::fem::FEMMetaData::NODE_RANK + 0: return Vertex;
      case stk::mesh::fem::FEMMetaData::NODE_RANK + 1: return Edge;
      case stk::mesh::fem::FEMMetaData::NODE_RANK + 2: return Face;
      case stk::mesh::fem::FEMMetaData::NODE_RANK + 3: return Volume;
    }
    return Vertex;
}
}

typedef stk::mesh::Field<double,stk::mesh::Cartesian>            CartesianField ;


/********************************************************
* Constructors                                          *
********************************************************/
STKMeshElement::STKMeshElement() :
    d_dim   ( 0),
    d_rank  ( 0),
    d_mesh   ( 0),
    d_meshID ( 0),
    ptr_element ( 0)
{
    typeID  = STKMeshElementTypeID();
    element = 0;
    d_globalID = MeshElementID();
}
STKMeshElement::STKMeshElement(int dim, stk::mesh::Entity* STKmesh_element, unsigned int rank, MeshID meshID, const STKMesh* mesh) :
    d_dim   ( dim),
    d_rank  ( rank),
    d_mesh   ( mesh),
    d_meshID ( meshID),
    ptr_element ( STKmesh_element)
{
    typeID  = STKMeshElementTypeID();
    element = NULL;
    AMP_ASSERT(STKmesh_element!=NULL);
    unsigned int local_id   = ptr_element->identifier();
    unsigned int owner_rank = ptr_element->owner_rank();
    GeomType type           = geom_type(ptr_element->entity_rank());
    const bool is_local     = owner_rank==d_rank;
    d_globalID = MeshElementID(is_local, type, local_id, owner_rank, meshID);
}
STKMeshElement::STKMeshElement(int dim, AMP::shared_ptr< stk::mesh::Entity > STKmesh_element,
                               unsigned int rank, MeshID meshID, const STKMesh* mesh) :
    d_dim   ( dim),
    d_rank  ( rank),
    d_mesh   ( mesh),
    d_meshID ( meshID),
    ptr_element ( STKmesh_element.get())
{
    typeID  = STKMeshElementTypeID();
    element = NULL;
    AMP_ASSERT(STKmesh_element.get()!=NULL);
    unsigned int local_id   = ptr_element->identifier();
    unsigned int owner_rank = ptr_element->owner_rank();
    GeomType type           = geom_type(ptr_element->entity_rank());
    const bool is_local     = owner_rank==d_rank;
    d_globalID = MeshElementID(is_local,type,local_id,owner_rank,meshID);
}
STKMeshElement::STKMeshElement(const STKMeshElement& rhs) :
    d_dim   ( rhs.d_dim),
    d_rank  ( rhs.d_rank),
    d_mesh   ( rhs.d_mesh),
    d_meshID ( rhs.d_meshID),
    ptr_element ( rhs.ptr_element)
{
    typeID  = STKMeshElementTypeID();
    element = rhs.element;
    d_globalID = rhs.d_globalID;
}

STKMeshElement& STKMeshElement::operator=(const STKMeshElement& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    this->typeID      = STKMeshElementTypeID();
    this->element     = 0;
    this->d_globalID  = rhs.d_globalID;
    this->d_dim       = rhs.d_dim;
    this->ptr_element = rhs.ptr_element;
    this->d_mesh      = rhs.d_mesh;
    this->d_meshID    = rhs.d_meshID;
    return *this;
}


/****************************************************************
* De-constructor                                                *
****************************************************************/
STKMeshElement::~STKMeshElement()
{
    ptr_element = NULL;
}


/****************************************************************
* Function to clone the element                                 *
****************************************************************/
MeshElement* STKMeshElement::clone() const
{
    return new STKMeshElement(*this);
}


/****************************************************************
* Function to get the elements composing the current element    *
****************************************************************/
std::vector<MeshElement> STKMeshElement::getElements(const GeomType type) const
{
    AMP_INSIST(type<=d_globalID.type(),"sub-elements must be of a smaller or equivalent type");
    std::vector<MeshElement> children(0);
    stk::mesh::Entity* elem = ptr_element;
    if ( d_globalID.type()==Vertex ) {
        // A vertex does not have children, return itself
        if ( type!=Vertex ) AMP_ERROR("A vertex is the base element and cannot have and sub-elements");
        children.resize(1);
        children[0] = *this;
    } else if ( type==d_globalID.type() ) {
        // Return the children of the current element
        const stk::mesh::PairIterRelation r = elem->relations(elem->entity_rank());
        if ( r.first != r.second ) {
            children.resize(r.second-r.first);
            unsigned c=0;
            for (stk::mesh::PairIterRelation::iterator i=r.first; i!=r.second; ++i,++c)
                children[c] = STKMeshElement( d_dim, i->entity(), d_rank, d_meshID, d_mesh );
        } else {
            children.resize(1);
            children[0] = *this;
        }
    } else if ( type==Vertex ) {
        // Return the nodes of the current element
        const stk::mesh::PairIterRelation r = elem->node_relations();
        children.resize(r.second-r.first);
        int c=0;
        for (stk::mesh::PairIterRelation::iterator i=r.first; i!=r.second; ++i,++c)
            children[c] = STKMeshElement( d_dim, i->entity(), d_rank, d_meshID, d_mesh );
    } else {
        // Return the children
        stk::mesh::PairIterRelation r;
        if ( type==Edge )
            r = elem->relations(d_mesh->d_STKMeshMeta->edge_rank());
        else if ( type==Face )
            r = elem->relations(d_mesh->d_STKMeshMeta->face_rank());
        else
            AMP_ERROR("Internal error");
        children.resize(r.second-r.first);
        unsigned c = 0;
        for (stk::mesh::PairIterRelation::iterator i=r.first; i!=r.second; ++i,++c) {
            // We need to build a valid element
            stk::mesh::Entity * elem = i->entity();
            // Create the STKMeshElement
            children[c] = STKMeshElement( d_dim, elem, d_rank, d_meshID, d_mesh );
        }
    }
    return children;
}


/****************************************************************
* Function to get the neighboring elements                      *
****************************************************************/
std::vector< MeshElement::shared_ptr > STKMeshElement::getNeighbors() const
{
    std::vector< MeshElement::shared_ptr > neighbors(0);
    if ( d_globalID.type()==Vertex ) {
        // Return the neighbors of the current node
        std::vector< stk::mesh::Entity* > neighbor_nodes = d_mesh->getNeighborNodes( d_globalID );
        neighbors.resize(neighbor_nodes.size(), MeshElement::shared_ptr());
        for (size_t i=0; i<neighbor_nodes.size(); i++) {
            // There are no NULL neighbors
            AMP::shared_ptr<STKMeshElement> neighbor(new STKMeshElement( d_dim, neighbor_nodes[i], d_rank, d_meshID, d_mesh ) );
            neighbors[i] = neighbor;
        }
    } else if ( (int) d_globalID.type() == d_dim ) {
        // Return the neighbors of the current element
        stk::mesh::Entity* elem = ptr_element;

        stk::mesh::EntityVector adjacent_entities;
        const CellTopologyData* celltopology = stk::mesh::fem::get_cell_topology(*elem).getCellTopologyData();
        AMP_INSIST(celltopology,"No topology for element, can not find neighbors.");
        const stk::mesh::EntityRank subcell_rank = elem->entity_rank()-1;
        const unsigned   num_sides = celltopology->subcell_count[subcell_rank];
        for (unsigned id = 0; id < num_sides; ++id) {
            stk::mesh::EntityVector subcell_nodes;
            stk::mesh::EntityVector adj_entities;
            const CellTopologyData * subcell_topology = stk::mesh::fem::get_subcell_nodes(*elem,
                                                                                          subcell_rank,
                                                                                          id,
                                                                                          subcell_nodes);
            NULL_USE(subcell_topology);
            get_entities_through_relations(subcell_nodes, elem->entity_rank(), adj_entities);
            adjacent_entities.insert(adjacent_entities.end(), adj_entities.begin(), adj_entities.end());
        }
        std::sort(adjacent_entities.begin(), adjacent_entities.end());
        adjacent_entities.resize(std::unique(adjacent_entities.begin(), adjacent_entities.end()) - adjacent_entities.begin());
        adjacent_entities.erase (std::find  (adjacent_entities.begin(), adjacent_entities.end(), elem));

        neighbors.resize(adjacent_entities.size());
        for (size_t i=0; i<neighbors.size(); i++) {
            stk::mesh::Entity *neighbor_elem = adjacent_entities[i];
            AMP::shared_ptr<STKMeshElement> neighbor(new STKMeshElement( d_dim, neighbor_elem, d_rank, d_meshID, d_mesh ) );
            neighbors[i] = neighbor;
        }
    } else {
        // We constructed a temporary STKmesh object and do not have access to the neighbor info
    }
    return neighbors;
}


/****************************************************************
* Functions to get basic element properties                     *
****************************************************************/
double STKMeshElement::volume() const
{
    const unsigned numCells = 1;
    if ( d_globalID.type() == Vertex )
        AMP_ERROR("STKMeshElement::volume:volume is is not defined on Nodes");
    if ( d_globalID.type() == Edge ) {
        AMP_WARNING("STKMeshElement::volume:volume is is not defined on Edges");
        return 1;
    }
    if ( d_globalID.type() == Face ) {
        AMP_WARNING("STKMeshElement::volume:volume is is not defined on Faces");
        return 1;
    }
    const stk::mesh::Entity* elem = (stk::mesh::Entity*) ptr_element;

    const CartesianField *coordinates = d_mesh->d_STKMeshMeta->get_field<CartesianField>("coordinates");
    const stk::mesh::PairIterRelation elem_nodes = elem->node_relations();
    const unsigned numNodes = elem_nodes.end() - elem_nodes.begin();

    const stk::mesh::fem::CellTopology cell_topo = stk::mesh::fem::get_cell_topology(*elem); 

    Teuchos::RCP<Intrepid::Cubature<double> > myCub;
    try {
        const unsigned cubDegree = 4; 
        Intrepid::DefaultCubatureFactory<double> cubFactory;
        myCub = cubFactory.create(cell_topo, cubDegree);
    }
    catch (...) {
        AMP_ERROR("STKMeshElement::volume mesh contains elements that Intrepid doesn't support for quadrature.");
    }

    const unsigned numCubPoints = myCub->getNumPoints(); 
    Intrepid::FieldContainer<double> volume         (numCells);
    Intrepid::FieldContainer<double> onesLeft       (numCells, numCubPoints);
    Intrepid::FieldContainer<double> weightedMeasure(numCells, numCubPoints);
    Intrepid::FieldContainer<double> jacobian_det   (numCells, numCubPoints);
    Intrepid::FieldContainer<double> jacobian       (numCells, numCubPoints, d_dim, d_dim);

    Intrepid::FieldContainer<double> cub_points               (numCubPoints, d_dim);
    Intrepid::FieldContainer<double> cub_weights              (numCubPoints);
    Intrepid::FieldContainer<double> cellNodes      (numCells, numNodes, d_dim);

    for (unsigned i=0; i!=numNodes; ++i) {
        const stk::mesh::Entity *node = elem_nodes[i].entity();
        const double *X =  (const double*)stk::mesh::field_data(*coordinates, *node);
        for (int j=0; j<d_dim; j++) 
            cellNodes(0,i,j) = X[j];
    }

    myCub->getCubature(cub_points, cub_weights); 
    Intrepid::CellTools<double>::setJacobian   (jacobian,     cub_points, cellNodes, cell_topo);
    Intrepid::CellTools<double>::setJacobianDet(jacobian_det, jacobian);

    onesLeft.initialize(1.0);
    Intrepid::FunctionSpaceTools::computeCellMeasure<double>(weightedMeasure, jacobian_det, cub_weights);
    Intrepid::FunctionSpaceTools::integrate<double>(volume, onesLeft, weightedMeasure,  Intrepid::COMP_BLAS);
    const double v = volume(0);
    return v;
}
std::vector<double> STKMeshElement::coord() const
{
    if ( d_globalID.type() != Vertex )
        AMP_ERROR("coord is only defined for Nodes");
    stk::mesh::Entity* node = (stk::mesh::Entity*) ptr_element;

    CartesianField *coordinates = d_mesh->d_STKMeshMeta->get_field<CartesianField>("coordinates");
    double *X =  (double*)stk::mesh::field_data(*coordinates, *node);

    const std::vector<double> x(X, &X[d_dim]);
    return x;
}
std::vector<double> STKMeshElement::centroid() const
{
    if ( d_globalID.type()==Vertex ) return coord();
    stk::mesh::Entity* elem = (stk::mesh::Entity*) ptr_element;

    std::vector<double> x(d_dim,0.0);
    CartesianField *coordinates = d_mesh->d_STKMeshMeta->get_field<CartesianField>("coordinates");
    stk::mesh::PairIterRelation elem_nodes = elem->node_relations();
    for (stk::mesh::PairIterRelation::iterator i=elem_nodes.begin(); i!=elem_nodes.end(); ++i) {
        const stk::mesh::Entity *node = i->entity();
        const double *X =  (double*)stk::mesh::field_data(*coordinates, *node);
        for (int j=0; j<d_dim; j++) x[j] += X[j];
    }
    const unsigned len = elem_nodes.end()-elem_nodes.begin(); 
    for (int i=0; i<d_dim; i++) x[i] /= len;
    return x;
}
bool STKMeshElement::containsPoint( const std::vector<double> &pos, double TOL ) const
{
    AMP_ERROR("STKMeshElement::containPoint is is not defined");
    if ( d_globalID.type()==Vertex ) {
        //double dist = 0.0;
        std::vector<double> point = this->coord();
        double dist2 = 0.0;
        for (size_t i=0; i<point.size(); i++)
            dist2 += (point[i]-pos[i])*(point[i]-pos[i]);
        return dist2<=TOL*TOL;
    }
    //stk::mesh::Entity* elem = (stk::mesh::Entity*) ptr_element;
    //std::vector<double> point ;//(pos[0],pos[1],pos[2]);
    return false; //elem->contains_point(point,TOL);
}
bool STKMeshElement::isOnSurface() const
{
    GeomType type = d_globalID.type();
    MeshElement search = MeshElement(*this);
    if ( d_globalID.is_local() ) {
        std::vector<MeshElement> &data = *(d_mesh->d_localSurfaceElements[type]);
        size_t index = AMP::Utilities::findfirst( data, search );
        if ( index < data.size() ) {
            if ( d_mesh->d_localSurfaceElements[type]->operator[](index).globalID() == d_globalID )
                return true;
        }
    } else {
        std::vector<MeshElement> &data = *(d_mesh->d_ghostSurfaceElements[type]);
        size_t index = AMP::Utilities::findfirst( data, search );
        if ( index < data.size() ) {
            if ( d_mesh->d_ghostSurfaceElements[type]->operator[](index).globalID() == d_globalID )
                return true;
        }
    }
    return false;
}
bool STKMeshElement::isOnBoundary(int) const
{
    return isOnSurface();
}
bool STKMeshElement::isInBlock(int id) const
{
   AMP_ERROR("STKMeshElement::isInBlock is is not defined");
   GeomType type = d_globalID.type();
   bool in_block = false;
    if ( type==Vertex ) {
        // Entity is a libmesh node
        AMP_ERROR("isInBlock is not currently implimented for anything but elements");
    } else if ( (int)type==d_dim ) {
        // Entity is a libmesh node
        //stk::mesh::Entity* elem = (stk::mesh::Entity*) ptr_element;
        in_block = false;//elem->subdomain_id() == id;
    } else  {
        // All other entities are on the boundary iff all of their verticies are on the surface
        AMP_ERROR("isInBlock is not currently implimented for anything but elements");
    }
    return in_block;
}



} // Mesh namespace
} // AMP namespace
