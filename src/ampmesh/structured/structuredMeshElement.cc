#include "ampmesh/structured/structuredMeshElement.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Mesh {


// Create a unique id for this class
static unsigned int structuredMeshElementTypeID = TYPE_HASH(structuredMeshElement);

// Function to evaluate the magnitude of a cross product in 3d
double cross3magnitude( double a[3], double b[3] )
{
    double v[3];
    v[0] = a[1]*b[2] - a[2]*b[1];
    v[1] = a[2]*b[0] - a[0]*b[2];
    v[2] = a[0]*b[1] - a[1]*b[0];
    return std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}
// Function to evaluate the dot produce of a vector anda cross product in 3d ( a . ( b X c ) )
double dot3cross( double a[3], double b[3], double c[3] )
{
    double v[3];
    v[0] = b[1]*c[2] - b[2]*c[1];
    v[1] = b[2]*c[0] - b[0]*c[2];
    v[2] = b[0]*c[1] - b[1]*c[0];
    return a[0]*v[0] + a[1]*v[1] + a[2]*v[2];
}


/********************************************************
* Constructors                                          *
********************************************************/
structuredMeshElement::structuredMeshElement()
{
    typeID = structuredMeshElementTypeID;
    element = NULL;
    d_dim = 0;
    d_index = BoxMesh::MeshElementIndex();
    d_globalID = MeshElementID();
}
structuredMeshElement::structuredMeshElement(  BoxMesh::MeshElementIndex index, const AMP::Mesh::BoxMesh* mesh )
{
    typeID = structuredMeshElementTypeID;
    d_mesh = mesh;
    d_dim = d_mesh->getDim();
    AMP_ASSERT(d_dim>0&&d_dim<=3);
    d_index = index;
    unsigned int owner_rank;
    int myBoxSize[3]  = {1,1,1};
    int myBoxRange[6] = {0,0,0,0,0,0};
    d_mesh->getOwnerBlock( d_index, owner_rank, myBoxRange );
    for (int d=0; d<d_mesh->PhysicalDim; d++) {
        AMP_ASSERT(index.index[d]>=myBoxRange[2*d+0]&&index.index[d]<myBoxRange[2*d+1]);
        myBoxSize[d] = myBoxRange[2*d+1] - myBoxRange[2*d+0];
    }
    unsigned int local_id = (index.index[0]-myBoxRange[0]) + (index.index[1]-myBoxRange[2])*myBoxSize[0] +
        (index.index[2]-myBoxRange[4])*myBoxSize[0]*myBoxSize[1] + index.side*myBoxSize[0]*myBoxSize[1]*myBoxSize[2];
    bool is_local = (int)owner_rank == d_mesh->d_comm.getRank();
    d_globalID = MeshElementID(is_local,(GeomType)index.type,local_id,owner_rank,d_mesh->d_meshID);
}
structuredMeshElement::structuredMeshElement(const structuredMeshElement& rhs)
{
    typeID = structuredMeshElementTypeID;
    element = NULL;
    d_globalID = rhs.d_globalID;
    d_dim = rhs.d_dim;
    d_index = rhs.d_index;
    d_mesh = rhs.d_mesh;
}
structuredMeshElement& structuredMeshElement::operator=(const structuredMeshElement& rhs)
{
    if (this == &rhs) // protect against invalid self-assignment
        return *this;
    this->typeID = structuredMeshElementTypeID;
    this->element = NULL;
    this->d_globalID = rhs.d_globalID;
    this->d_dim = rhs.d_dim;
    this->d_index = rhs.d_index;
    this->d_mesh = rhs.d_mesh;
    return *this;
}


/****************************************************************
* De-constructor                                                *
****************************************************************/
structuredMeshElement::~structuredMeshElement()
{
}


/****************************************************************
* Function to clone the element                                 *
****************************************************************/
MeshElement* structuredMeshElement::clone() const
{
    return new structuredMeshElement(*this);
}


/****************************************************************
* Function to get the elements composing the current element    *
* We use a Canonical numbering system                           *
****************************************************************/
std::vector<MeshElement> structuredMeshElement::getElements(const GeomType type) const
{
    AMP_ASSERT(type<=d_dim);
    if ( type==d_globalID.type() )
        return std::vector<MeshElement>(1,MeshElement(*this));
    std::vector<BoxMesh::MeshElementIndex> index;
    if ( type==Vertex ) {
        // We want to get the verticies composing the elements
        if ( d_globalID.type()==(GeomType)d_dim ) {
            // We are dealing with a entity of type dim and want the verticies
            if ( d_dim==1 ) {
                index.resize(2);
                index[0] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]   );
                index[1] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1 );
            } else if ( d_dim==2 ) {
                index.resize(4);
                index[0] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1]   );
                index[1] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1]   );
                index[2] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1]+1 );
                index[3] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1]+1 );
            } else if ( d_dim==3 ) {
                index.resize(8);
                index[0] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1],   d_index.index[2]   );
                index[1] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1],   d_index.index[2]   );
                index[2] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1]+1, d_index.index[2]   );
                index[3] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1]+1, d_index.index[2]   );
                index[4] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1],   d_index.index[2]+1 );
                index[5] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1],   d_index.index[2]+1 );
                index[6] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1]+1, d_index.index[2]+1 );
                index[7] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1]+1, d_index.index[2]+1 );
            } else {
                AMP_ERROR("Dimension not supported yet");
            }
        } else if ( d_globalID.type()==Edge ) {
            index.resize(2,d_index);
            index[0].type = 0;
            index[1].type = 0;
            index[0].side = 0;
            index[1].side = 0;
            index[0].index[d_index.side] = d_index.index[d_index.side];
            index[1].index[d_index.side] = d_index.index[d_index.side]+1;
        } else if ( d_globalID.type()==Face ) {
            index.resize(4,d_index);
            if ( d_index.side==0 ) {
                index[0] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0], d_index.index[1],   d_index.index[2]   );
                index[1] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0], d_index.index[1]+1, d_index.index[2]   );
                index[2] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0], d_index.index[1]+1, d_index.index[2]+1 );
                index[3] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0], d_index.index[1],   d_index.index[2]+1 );
            } else if ( d_index.side==1 ) {
                index[0] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1], d_index.index[2]   );
                index[1] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1], d_index.index[2]   );
                index[2] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1], d_index.index[2]+1 );
                index[3] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1], d_index.index[2]+1 );
            } else if ( d_index.side==2 ) {
                index[0] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1],   d_index.index[2] );
                index[1] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1],   d_index.index[2] );
                index[2] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1]+1, d_index.index[2] );
                index[3] = BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1]+1, d_index.index[2] );
            } else {
                AMP_ERROR("Internal error");
            }
        } else if ( d_globalID.type()==Volume ) {
            AMP_ERROR("Not ready for dimensions > 3");
        } else {
            AMP_ERROR("Not finsihed");
        }
    } else if ( type==Edge ) {
        if ( d_globalID.type()==Face ) {
            int side = d_index.side;
            if ( d_dim==2 )
                side = 2;
            index.resize(4,d_index);
            if ( d_index.side==0 ) {
                // We are dealing with an x-face
                index[0] = BoxMesh::MeshElementIndex( Edge, 1, d_index.index[0], d_index.index[1],   d_index.index[2]   );
                index[1] = BoxMesh::MeshElementIndex( Edge, 2, d_index.index[0], d_index.index[1]+1, d_index.index[2]   );
                index[2] = BoxMesh::MeshElementIndex( Edge, 1, d_index.index[0], d_index.index[1],   d_index.index[2]+1 );
                index[3] = BoxMesh::MeshElementIndex( Edge, 2, d_index.index[0], d_index.index[1],   d_index.index[2]   );
            } else if ( d_index.side==1 ) {
                // We are dealing with an y-face
                index[0] = BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0],   d_index.index[1], d_index.index[2]   );
                index[1] = BoxMesh::MeshElementIndex( Edge, 2, d_index.index[0]+1, d_index.index[1], d_index.index[2]   );
                index[2] = BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0],   d_index.index[1], d_index.index[2]+1 );
                index[3] = BoxMesh::MeshElementIndex( Edge, 2, d_index.index[0],   d_index.index[1], d_index.index[2]   );
            } else if ( d_index.side==2 ) {
                // We are dealing with an z-face
                index[0] = BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0],   d_index.index[1],   d_index.index[2] );
                index[1] = BoxMesh::MeshElementIndex( Edge, 1, d_index.index[0]+1, d_index.index[1],   d_index.index[2] );
                index[2] = BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0],   d_index.index[1]+1, d_index.index[2] );
                index[3] = BoxMesh::MeshElementIndex( Edge, 1, d_index.index[0],   d_index.index[1],   d_index.index[2] );
            } else {
                AMP_ERROR("Internal error");
            }
        } else if ( d_globalID.type()==Volume ) {
            AMP_ASSERT(d_index.side==0);
            index.resize(12,d_index);
            index[0] = BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0],   d_index.index[1],   d_index.index[2]   );
            index[1] = BoxMesh::MeshElementIndex( Edge, 1, d_index.index[0]+1, d_index.index[1],   d_index.index[2]   );
            index[2] = BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0],   d_index.index[1]+1, d_index.index[2]   );
            index[3] = BoxMesh::MeshElementIndex( Edge, 1, d_index.index[0],   d_index.index[1],   d_index.index[2]   );
            index[4] = BoxMesh::MeshElementIndex( Edge, 2, d_index.index[0],   d_index.index[1],   d_index.index[2]   );
            index[5] = BoxMesh::MeshElementIndex( Edge, 2, d_index.index[0]+1, d_index.index[1],   d_index.index[2]   );
            index[6] = BoxMesh::MeshElementIndex( Edge, 2, d_index.index[0]+1, d_index.index[1]+1, d_index.index[2]   );
            index[7] = BoxMesh::MeshElementIndex( Edge, 2, d_index.index[0],   d_index.index[1],   d_index.index[2]   );
            index[8] = BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0],   d_index.index[1],   d_index.index[2]+1 );
            index[9] = BoxMesh::MeshElementIndex( Edge, 1, d_index.index[0]+1, d_index.index[1],   d_index.index[2]+1 );
            index[10] = BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0],   d_index.index[1]+1, d_index.index[2]+1 );
            index[11] = BoxMesh::MeshElementIndex( Edge, 1, d_index.index[0],   d_index.index[1],   d_index.index[2]+1 );
        } else {
            AMP_ERROR("Dimensions > 3 are not supported yet");
        }
    } else if ( type==Face ) {
        if ( d_globalID.type()==Volume ) {
            index.resize(6,d_index);
            index[0] = BoxMesh::MeshElementIndex( Face, 1, d_index.index[0],   d_index.index[1],   d_index.index[2]   );
            index[1] = BoxMesh::MeshElementIndex( Face, 0, d_index.index[0]+1, d_index.index[1],   d_index.index[2]   );
            index[2] = BoxMesh::MeshElementIndex( Face, 1, d_index.index[0],   d_index.index[1]+1, d_index.index[2]   );
            index[3] = BoxMesh::MeshElementIndex( Face, 0, d_index.index[0],   d_index.index[1],   d_index.index[2]   );
            index[4] = BoxMesh::MeshElementIndex( Face, 2, d_index.index[0],   d_index.index[1],   d_index.index[2]   );
            index[5] = BoxMesh::MeshElementIndex( Face, 2, d_index.index[0],   d_index.index[1],   d_index.index[2]+1 );
        } else {
            AMP_ERROR("Dimensions > 3 are not supported yet");
        }
    } else if ( type==Volume ) {
        AMP_ERROR("Dimensions > 3 are not supported yet");
    } else  {
        AMP_ERROR("Not finished");
    }
    // Fix any elements that are beyond a periodic boundary
    for (int d=0; d<d_dim; d++) {
        if ( d_mesh->d_isPeriodic[d] ) {
            int size = d_mesh->d_size[d];;
            for (size_t i=0; i<index.size(); i++) {
                if ( index[i].index[d]<0 )
                    index[i].index[d] += size;
                else if ( index[i].index[d]>=size )
                    index[i].index[d] -= size;
            }
        }
    }    
    // Get the elements
    std::vector<MeshElement> elements(index.size());
    for (size_t i=0; i<index.size(); i++)
        elements[i] = structuredMeshElement( index[i], d_mesh );
    return elements;
}


/****************************************************************
* Function to get the neighboring elements                      *
****************************************************************/
std::vector<MeshElement::shared_ptr> structuredMeshElement::getNeighbors() const
{
    std::vector<BoxMesh::MeshElementIndex> index;
    if ( d_globalID.type()==Vertex ) {
        // Get the list of neighbor nodex (there are no null neighbors)
        // The node neighbors are the list of nodes that share any element
        index.reserve(27);
        if ( d_dim==1 ) {
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]-1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1 ) );
        } else if ( d_dim==2 ) {
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]-1, d_index.index[1]-1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1]-1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1]-1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]-1, d_index.index[1]   ) );
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1]   ) );
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]-1, d_index.index[1]+1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0],   d_index.index[1]+1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+1, d_index.index[1]+1 ) );
        } else if ( d_dim==3 ) {
            for (int k=-1; k<=1; k++) {
                for (int j=-1; j<=1; j++) {
                    for (int i=-1; i<=1; i++) {
                        if ( i==0 && j==0 && k==0 )
                            continue;
                        index.push_back( BoxMesh::MeshElementIndex( Vertex, 0, d_index.index[0]+i, d_index.index[1]+j, d_index.index[2]+k ) );
                    }
                }
            }
        } else {
            AMP_ERROR("Dimension not supported yet");
        }
    } else if ( d_globalID.type()==Edge ) {
        if ( d_dim==1 ) {
            index.push_back( BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0]-1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Edge, 0, d_index.index[0]+1 ) );
        } else {
            // Edge neighbors in dimensions > 1 are not supported yet
        }
    } else if ( d_globalID.type()==Face ) {
        if ( d_dim==2 ) {
            index.push_back( BoxMesh::MeshElementIndex( Face, 0, d_index.index[0],   d_index.index[1]-1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Face, 0, d_index.index[0]+1, d_index.index[0]   ) );
            index.push_back( BoxMesh::MeshElementIndex( Face, 0, d_index.index[0],   d_index.index[1]+1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Face, 0, d_index.index[0]-1, d_index.index[1]   ) );
        } else {
            // Face neighbors in dimensions > 2 are not supported yet
        }
    } else if ( d_globalID.type()==Volume ) {
        if ( d_dim==3 ) {
            index.push_back( BoxMesh::MeshElementIndex( Volume, 0, d_index.index[0],   d_index.index[1]-1, d_index.index[2]   ) );
            index.push_back( BoxMesh::MeshElementIndex( Volume, 0, d_index.index[0]+1, d_index.index[1],   d_index.index[2]   ) );
            index.push_back( BoxMesh::MeshElementIndex( Volume, 0, d_index.index[0],   d_index.index[1]+1, d_index.index[2]   ) );
            index.push_back( BoxMesh::MeshElementIndex( Volume, 0, d_index.index[0]-1, d_index.index[1],   d_index.index[2]   ) );
            index.push_back( BoxMesh::MeshElementIndex( Volume, 0, d_index.index[0],   d_index.index[1],   d_index.index[2]-1 ) );
            index.push_back( BoxMesh::MeshElementIndex( Volume, 0, d_index.index[0],   d_index.index[1],   d_index.index[2]+1 ) );
        } else {
            // Volume neighbors in dimensions > 3 are not supported yet
        }
    } else {
        AMP_ERROR("Unknown entity type");
    }
    // Get the neighbor elements
    std::vector<MeshElement::shared_ptr> neighbors;
    neighbors.reserve(index.size());
    bool periodic[3];
    int size[3];
    for (int d=0; d<d_dim; d++) {
        periodic[d] = d_mesh->d_isPeriodic[d];
        size[d] = d_mesh->d_size[d];
    }
    for (size_t i=0; i<index.size(); i++) {
        bool in_mesh = true;
        for (int d=0; d<d_dim; d++) {
            if ( periodic[d] ) {
                if ( index[i].index[d]<0 )
                    index[i].index[d] += size[d];
                if ( index[i].index[d]>=size[d] )
                    index[i].index[d] -= size[d];
            } else {
                if ( index[i].index[d]<0 )
                    in_mesh = false;
                if ( d_globalID.type()==d_dim ) {
                    if ( index[i].index[d]>=size[d]  )
                        in_mesh = false;
                } else {
                    if ( index[i].index[d]>size[d]  )
                        in_mesh = false;
                }
            }
        }
        if ( in_mesh )
            neighbors.push_back( MeshElement::shared_ptr( new structuredMeshElement( index[i], d_mesh ) ) );
        else if ( d_globalID.type()!=Vertex ) 
            neighbors.push_back( MeshElement::shared_ptr() );
    }
    return neighbors;
}


/****************************************************************
* Functions to get the element volume                           *
****************************************************************/
double structuredMeshElement::volume() const
{
    if ( d_globalID.type() == Vertex ) {
        AMP_ERROR("volume is is not defined Nodes");
    }
    std::vector<MeshElement> nodes = this->getElements(Vertex);
    if ( d_globalID.type() == Edge ) {
        AMP_ASSERT(nodes.size()==2);
        std::vector<double> x1 = nodes[0].coord();
        std::vector<double> x2 = nodes[1].coord();
        double dist2 = 0.0;
        for (int i=0; i<d_dim; i++)
            dist2 += (x1[i]-x2[i])*(x1[i]-x2[i]);
        return sqrt(dist2);
    } else if ( d_globalID.type() == Face ) {
        // Use 2x2 quadrature to approximate the surface area. See for example,
        // Y. Zhang, C. Bajaj, G. Xu. Surface Smoothing and Quality Improvement
        // of Quadrilateral/Hexahedral Meshes with Geometric Flow. The special 
        // issue of the Journal Communications in Numerical Methods in
        // Engineering (CNME), submitted as an invited paper, 2006.
        // http://www.ices.utexas.edu/~jessica/paper/quadhexgf/quadhex_geomflow_CNM.pdf
        std::vector<double> x1 = nodes[0].coord();
        std::vector<double> x2 = nodes[1].coord();
        std::vector<double> x3 = nodes[2].coord();
        std::vector<double> x4 = nodes[3].coord();
        double AB[3]={0,0,0}, AC[3]={0,0,0}, AD[3]={0,0,0}, AC_AB_AD[3]={0,0,0};
        for (int i=0; i<d_dim; i++) {
            AC[i] = x3[i]-x1[i];                // Vector pointing from A to C
            AB[i] = x2[i]-x1[i];                // Vector pointing from A to B
            AD[i] = x4[i]-x1[i];                // Vector pointing from A to D
            AC_AB_AD[i] = AC[i]-AB[i]-AD[i];    // The diagonal vector minus the side vectors
        }
        if ( AC_AB_AD[0]==0 && AC_AB_AD[1]==0 && AC_AB_AD[2]==0 ) {
            // The points are co-planar
            return cross3magnitude(AB,AD);
        } else {
            const double q[2] = { 0.5-std::sqrt(3.0)/6.0, 0.5 + std::sqrt(3.0)/6.0 };
            double vol = 0.0;
            double v1[3], v2[3];
            for (unsigned int i=0; i<2; ++i) {
            	for (unsigned int j=0; j<2; ++j) {
                    v1[0] = AB[0] + q[i]*AC_AB_AD[0];
                    v1[1] = AB[1] + q[i]*AC_AB_AD[1];
                    v1[2] = AB[2] + q[i]*AC_AB_AD[2];
                    v2[0] = AD[0] + q[j]*AC_AB_AD[0];
                    v2[1] = AD[1] + q[j]*AC_AB_AD[1];
                    v2[2] = AD[2] + q[j]*AC_AB_AD[2];
                    vol += cross3magnitude(v1,v2);
                }
            }
            return 0.25*vol;
        }
    } else if ( d_globalID.type() == Volume ) {
        // Compute the volume of the tri-linear hex by splitting it
        // into 6 sub-pyramids and applying the formula in:
        // "Calculation of the Volume of a General Hexahedron
        // for Flow Predictions", AIAA Journal v.23, no.6, 1984, p.954-
        static const unsigned char sub_pyr[6][4] = {
              {0, 3, 2, 1}, {6, 7, 4, 5}, {0, 1, 5, 4},
              {3, 7, 6, 2}, {0, 4, 7, 3}, {1, 2, 6, 5} };
        // The centroid is a convenient point to use
        // for the apex of all the pyramids.
        std::vector<double> R = this->centroid();
        std::vector<std::vector<double> > points(8);
        for (int i=0; i<8; i++)
            points[i] = nodes[i].coord();
        int pyr_base[4];
        double vol=0.0;
        // Compute the volume using 6 sub-pyramids
        for (unsigned int n=0; n<6; ++n) {
            // Set the nodes of the pyramid base
            for (unsigned int i=0; i<4; ++i)
            	pyr_base[i] = sub_pyr[n][i];
            // Compute diff vectors
            double a[3], b[3], c[3], d[3], e[3];
            for (int i=0; i<3; i++) {
                a[i] = points[pyr_base[0]][i] - R[i];
                b[i] = points[pyr_base[1]][i] - points[pyr_base[3]][i];
                c[i] = points[pyr_base[2]][i] - points[pyr_base[0]][i];
                d[i] = points[pyr_base[3]][i] - points[pyr_base[0]][i];
                e[i] = points[pyr_base[1]][i] - points[pyr_base[0]][i];
            }
            // Compute pyramid volume
            double sub_vol = (1.0/6.0)*dot3cross(a,b,c) + (1.0/12.0)*dot3cross(c,d,e);
            AMP_ASSERT(sub_vol>0.0);
            vol += sub_vol;
        }
        return vol;
    }
    AMP_ERROR("Internal error");
    return 0.0;
}


/****************************************************************
* Function to get the node coordinates                          *
****************************************************************/
std::vector<double> structuredMeshElement::coord() const
{
    AMP_ASSERT(d_globalID.type()==Vertex);
    size_t pos = AMP::Utilities::findfirst(d_mesh->d_index,d_index);
    AMP_ASSERT(d_mesh->d_index[pos]==d_index);
    std::vector<double> coord((size_t)d_dim,0.0);
    for (int i=0; i<d_dim; i++)
        coord[i] = d_mesh->d_coord[i][pos];
    return coord;
}


/****************************************************************
* Misc functions                                                *
****************************************************************/
bool structuredMeshElement::containsPoint( const std::vector<double> &pos, double TOL ) const
{
    AMP_ERROR("Not finsihed");
    return false;
}
bool structuredMeshElement::isOnSurface() const
{
    bool on_surface = false;
    for (int d=0; d<d_dim; d++) {
        if ( d_mesh->d_isPeriodic[d] )
            continue;
        int size = (int) d_mesh->d_size[d];
        if ( d_globalID.type()==d_mesh->GeomDim ) {
            // We are dealing with the highest level geometric entity
            if ( d_index.index[d]==0 || d_index.index[d]==size-1 )
                on_surface = true;
        } else if ( d_globalID.type()==Vertex ) {
            // We are dealing with a vertex
            if ( d_index.index[d]==0 || d_index.index[d]==size )
                on_surface = true;
        } else if ( d_globalID.type()==Edge ) {
            // We are dealing with a vertex
            if ( ( d_index.index[d]==0 || d_index.index[d]==size ) && d_index.side!=d )
                on_surface = true;
        } else if ( d_globalID.type()==Face ) {
            // We are dealing with a vertex
            if ( ( d_index.index[d]==0 || d_index.index[d]==size ) && d_index.side==d )
                on_surface = true;
        } else {
            AMP_ERROR("Internal error (dim>3?)");
        }
    }
    return on_surface;
}
bool structuredMeshElement::isOnBoundary(int id) const
{
    std::map<std::pair<int,GeomType>,std::vector<BoxMesh::ElementIndexList> >::const_iterator iterator;
    iterator = d_mesh->d_id_list.find( std::pair<int,GeomType>( id, (GeomType)d_index.type ) );
    if ( iterator==d_mesh->d_id_list.end() )
        return false;
    for (size_t i=0; i<iterator->second.size(); i++) {
        if ( iterator->second[i]->size()==0 )
            continue;
        size_t j = AMP::Utilities::findfirst( *(iterator->second[i]), d_index );
        if ( j==iterator->second[i]->size() ) { j--; }
        if ( iterator->second[i]->operator[](j) == d_index )
            return true;
    }
    return false;
}
bool structuredMeshElement::isInBlock(int id) const
{
    if ( id==0 )
        return true;
    return false;
}


} // Mesh namespace
} // AMP namespace

