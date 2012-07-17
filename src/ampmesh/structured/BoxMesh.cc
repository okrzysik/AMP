#include "ampmesh/structured/BoxMesh.h"
#include "ampmesh/structured/structuredMeshElement.h"
#include "ampmesh/structured/structuredMeshIterator.h"
#include "ampmesh/MultiIterator.h"

#include "utils/Utilities.h"
#ifdef USE_AMP_VECTORS
    #include "vectors/Vector.h"
    #include "vectors/Variable.h"
    #include "vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
    #include "discretization/DOF_Manager.h"
    #include "discretization/simpleDOF_Manager.h"
#endif
#include "utils/ProfilerApp.h"

#include <iostream>
#include <algorithm>

namespace AMP {
namespace Mesh {


/****************************************************************
* Constructor                                                   *
****************************************************************/
BoxMesh::BoxMesh( const MeshParameters::shared_ptr &params_in ):
    Mesh(params_in)
{
    PROFILE_START("Constructor");
    for (int d=0; d<3; d++) {
        d_size[d] = 1;
        d_isPeriodic[d] = false;
        d_numBlocks[d] = 1;
    }
    // Check for valid inputs
    AMP_INSIST(params.get(),"Params must not be null");
    AMP_INSIST(d_comm!=AMP_MPI(AMP_COMM_NULL),"Communicator must be set");
    AMP_INSIST(d_db.get(),"Database must exist");
    // Get mandatory fields from the database
    AMP_INSIST(d_db->keyExists("dim"),"Field 'dim' must exist in database'");
    AMP_INSIST(d_db->keyExists("Generator"),"Field 'Generator' must exist in database'");
    AMP_INSIST(d_db->keyExists("Size"),"Field 'Size' must exist in database'");
    AMP_INSIST(d_db->keyExists("Range"),"Field 'Range' must exist in database'");
    PhysicalDim = d_db->getInteger("dim");
    GeomDim = (GeomType) PhysicalDim;
    std::string generator = d_db->getString("Generator");
    std::vector<int> size = d_db->getIntegerArray("Size");
    std::vector<double> range = d_db->getDoubleArray("Range");
    d_max_gcw = d_db->getIntegerWithDefault("GCW",2);
    for (size_t d=0; d<size.size(); d++)
        AMP_INSIST(size[d]>0,"All dimensions must have a size > 0");
    // Create the logical mesh
    AMP_ASSERT(PhysicalDim<=3);
    if ( generator.compare("cube")==0 ) {
        AMP_INSIST(size.size()==PhysicalDim,"Size of field 'Size' must match dim");
        for (int d=0; d<PhysicalDim; d++)
            d_size[d] = size[d];
    } else if ( generator.compare("cylinder")==0 ) {
        AMP_INSIST(PhysicalDim==3,"cylinder generator requires a 3d mesh");
        AMP_INSIST(size.size()==2,"Size of field 'Size' must be of size 2");
        d_size[0] = 2*size[0];
        d_size[1] = 2*size[0];
        d_size[2] = size[1];
    } else if ( generator.compare("tube")==0 ) {
        AMP_INSIST(PhysicalDim==3,"tube generator requires a 3d mesh");
        AMP_INSIST(size.size()==PhysicalDim,"Size of field 'Size' must match dim");
        for (int d=0; d<PhysicalDim; d++)
            d_size[d] = size[d];
        d_isPeriodic[1] = true;    // We will use the logical mesh (r,theta,z), so theta is periodic
    } else if ( generator.compare("sphere")==0 ) {
        AMP_INSIST(PhysicalDim==3,"sphere generator requires a 3d mesh");
        d_isPeriodic[1] = true;    // We will use the logical mesh (r,theta,phi), so theta and phi are periodic
        d_isPeriodic[2] = true;
        AMP_ERROR("Not finished");
    } 
    // Create the load balance
    for (int d=0; d<PhysicalDim; d++)
        d_numBlocks[d] = 1;
    if ( d_comm.getSize()==1 ) {
        // We are dealing with a serial mesh (do nothing to change the local box sizes)
    } else {
        // We are dealing with a parallel mesh
        // First, get the prime factors for number of processors and divide the dimensions
        std::vector<int> factors = AMP::Utilities::factor(d_comm.getSize());
        std::vector<int> div(PhysicalDim,1);
        while ( factors.size() > 0 ) {
            int d = -1;
            double v = -1;
            for (int i=0; i<PhysicalDim; i++) {
                double tmp = ((double)d_size[i])/((double)d_numBlocks[i]);
                if ( tmp > v ) {
                    d = i;
                    v = tmp;
                }
            }
            d_numBlocks[d] *= factors[factors.size()-1];
            factors.resize(factors.size()-1);
        }
    }
    AMP_INSIST(d_numBlocks[0]*d_numBlocks[1]*d_numBlocks[2]==d_comm.getSize(),
        "Error creating proper load balance");
    // Initialize the logical mesh
    initialize();
    // Create the appropriate mesh coordinates (and modify the id sets if necessary)
    if ( generator.compare("cube")==0 ) {
        AMP_INSIST(range.size()==2*PhysicalDim,"Range must be 2*dim for cube generator");
        fillCartesianNodes( PhysicalDim, &d_size[0], &range[0], d_index, d_coord );
    } else if ( generator.compare("tube")==0 ) {
        AMP_INSIST(range.size()==4,"Range must be 1x4 for tube generator");
        // Change the surface ids to match the standard ids
        // 0 - 8: Inner surface
        // 1 - 4: Outer surface
        // 4 - 2: Bottom surface
        // 5 - 1: Top surface
        AMP_ASSERT(d_ids.size()==4);
        d_ids[0] = 1;
        d_ids[1] = 2;
        d_ids[2] = 4;
        d_ids[3] = 8;
        std::map<std::pair<int,GeomType>,std::vector<ElementIndexList> >::const_iterator it;
        std::map<std::pair<int,GeomType>,std::vector<ElementIndexList> > new_ids;
        for (it=d_id_list.begin(); it!=d_id_list.end(); ++it) {
            std::pair<int,GeomType> id = it->first;
            std::vector<ElementIndexList> list = it->second;
            if ( id.first==0 )
                id.first = 8;
            else if ( id.first==1 )
                id.first = 4;
            else if ( id.first==4 )
                id.first = 2;
            else if ( id.first==5 )
                id.first = 1;
            else
                AMP_ERROR("Unexpected id");
            std::pair<std::pair<int,GeomType>,std::vector<ElementIndexList> > tmp( id, list );
            new_ids.insert( tmp );
        }
        d_id_list = new_ids;
        // Create the coordinates (currently the points lie in [0,1])
        double *x = &d_coord[0][0];
        double *y = &d_coord[1][0];
        double *z = &d_coord[2][0];
        const double pi = 3.141592653589793116;
        for (size_t i=0; i<d_coord[0].size(); i++) {
            double r = range[0] + x[i]*(range[1]-range[0]);
            double theta = 2*pi*y[i];
            x[i] = r*cos(theta);
            y[i] = r*sin(theta);
            z[i] = range[2] + z[i]*(range[3]-range[2]);
        }
    } else if ( generator.compare("cylinder")==0 ) {
        AMP_INSIST(range.size()==3,"Range must be 1x3 for cylinder generator");
        // Change the surface ids to match the standard ids
        // 0,1,2,3 - 4: Outer surface
        // 4 - 2: Bottom surface
        // 5 - 1: Top surface
        AMP_ASSERT(d_ids.size()==6);
        d_ids.resize(3);
        d_ids[0] = 1;
        d_ids[1] = 2;
        d_ids[2] = 4;
        std::map<std::pair<int,GeomType>,std::vector<ElementIndexList> >::const_iterator it1, it2;
        std::map<std::pair<int,GeomType>,std::vector<ElementIndexList> > new_ids;
        for (it1=d_id_list.begin(); it1!=d_id_list.end(); ++it1) {
            std::pair<int,GeomType> id = it1->first;
            std::vector<ElementIndexList> list = it1->second;
            if ( id.first==0 || id.first==1 || id.first==2 || id.first==3 )
                id.first = 4;
            else if ( id.first==4 )
                id.first = 2;
            else if ( id.first==5 )
                id.first = 1;
            else
                AMP_ERROR("Unexpected id");
            it2 = new_ids.find( id );
            if ( it2==new_ids.end() ) {
                std::pair<std::pair<int,GeomType>,std::vector<ElementIndexList> > tmp( id, list );
                new_ids.insert( tmp );
            } else {
                AMP_ASSERT(it2->second.size()==list.size());
                for (size_t i=0; i<list.size(); i++) {
                    it2->second[i]->reserve( it2->second[i]->size()+list[i]->size() );
                    for (size_t j=0; j<list[i]->size(); j++)
                        it2->second[i]->push_back( list[i]->operator[](j) );
                    AMP::Utilities::unique( *(it2->second[i]) );
                }
            }
        }
        d_id_list = new_ids;
        // Create the coordinates (currently the points lie in [0,1])
        double r1 = range[0];
        double *x = &d_coord[0][0];
        double *y = &d_coord[1][0];
        double *z = &d_coord[2][0];
        // Perform the mapping for the circle
        const double sqrt2 = 1.41421356237;
        for (size_t i=0; i<d_coord[0].size(); i++) {
            // map [0,1] x [0,1] to circle of radius r1
            double xc = 2*(x[i]-0.5);                   // Change domain to [-1,1]
            double yc = 2*(y[i]-0.5);                   // Change domain to [-1,1]
            if ( fabs(xc)<1e-12 && fabs(yc)<1e-12 ) {
                // We are dealing with the center point
                x[i] = 0.0;
                y[i] = 0.0;
                continue;
            }
            double d = std::max(fabs(xc),fabs(yc));     // value on diagonal of computational grid
            double D = r1*d/sqrt2;                      // mapping d to D(d)
            double R = r1;                              // mapping d to R(d) (alternative is R = r1*d)
            double center = D - sqrt(R*R-D*D);
            double xp = D/d*fabs(xc);
            double yp = D/d*fabs(yc);
            if ( fabs(yc)>=fabs(xc) )
                yp = center + sqrt(R*R-xp*xp);
            if ( fabs(xc)>=fabs(yc) )
                xp = center + sqrt(R*R-yp*yp);
            if ( xc<0.0 )
                xp = -xp;
            if ( yc<0.0 )
                yp = -yp;
            x[i] = xp;
            y[i] = yp;
        }
        // Perform the mapping for z
        for (size_t i=0; i<d_coord[0].size(); i++)
            z[i] = range[1] + z[i]*(range[2]-range[1]);
    } else { 
        AMP_ERROR("Unknown generator");
    }
    // Fill in the final info for the mesh
    AMP_INSIST(d_db->keyExists("MeshName"),"MeshName must exist in input database");
    d_name = d_db->getString("MeshName");
    d_box_local = std::vector<double>(2*PhysicalDim);
    for (int d=0; d<PhysicalDim; d++) {
        d_box_local[2*d+0] = 1e100;
        d_box_local[2*d+1] = -1e100;
        for (size_t i=0; i<d_coord[d].size(); i++) {
            if ( d_coord[d][i]<d_box_local[2*d+0])
                d_box_local[2*d+0] = d_coord[d][i];
            if ( d_coord[d][i]>d_box_local[2*d+1])
                d_box_local[2*d+1] = d_coord[d][i];
        }
    }
    d_box = std::vector<double>(PhysicalDim*2);
    for (int i=0; i<PhysicalDim; i++) {
        d_box[2*i+0] = d_comm.minReduce( d_box_local[2*i+0] );
        d_box[2*i+1] = d_comm.maxReduce( d_box_local[2*i+1] );
    } 
    // Displace the mesh
    std::vector<double> displacement(PhysicalDim,0.0);
    if ( d_db->keyExists("x_offset") )
        displacement[0] = d_db->getDouble("x_offset");
    if ( d_db->keyExists("y_offset") )
        displacement[1] = d_db->getDouble("y_offset");
    if ( d_db->keyExists("z_offset") )
        displacement[2] = d_db->getDouble("z_offset");
    bool test = false;
    for (size_t i=0; i<displacement.size(); i++) {
        if ( displacement[i] != 0.0 )
           test = true;
    }        
    if ( test )
        displaceMesh(displacement);
    PROFILE_STOP("Constructor");
}


/****************************************************************
* Initialize the mesh                                           *
****************************************************************/
void BoxMesh::initialize()
{
    PROFILE_START("initialize");
    // Compute the element indicies for all local and ghost elements
    for (int d=0; d<=PhysicalDim; d++) 
        d_elements[d] = std::vector<ElementIndexList>(d_max_gcw+1);
    std::vector<int> range = getLocalBlock(d_comm.getRank());
    // First get the list of owned elements of each type
    PROFILE_START("create_owned_elements");
    size_t N_localElements = 1;
    for (int d=0; d<PhysicalDim; d++) {
        int local_size = range[2*d+1]-range[2*d+0];
        AMP_ASSERT(local_size>0);
        N_localElements *= (size_t) local_size;
    }
    for (int d=0; d<=PhysicalDim; d++) {
        d_elements[d][0].reset( new std::vector<MeshElementIndex>() );
        if ( d==0 || d==PhysicalDim ) {
            d_elements[d][0]->reserve( N_localElements );
        } else if ( d==1 ) {
            d_elements[d][0]->reserve( 6*N_localElements );     
        } else if ( d==2 ) {
            d_elements[d][0]->reserve( 3*N_localElements );
        } else {
            AMP_ERROR("Internal error");
        }
    }
    for (int d=0; d<=PhysicalDim; d++) {    // Loop through the geometric entities
        int numSides = 0;
        if ( d==0 || d==PhysicalDim ) {
            numSides = 1;
        } else if ( PhysicalDim==2 ) {
            numSides = 2;
        } else if ( PhysicalDim==3 ) {
            numSides = 3;
        }
        AMP_ASSERT(numSides>0);
        for (int s=0; s<numSides; s++) {
            // Extend the range for the last box to include the physical boundary
            std::vector<int> range2 = range;
            if ( d==PhysicalDim ) {
                // We are dealing with an element, do nothing
            } else if ( d==0 ) {
                // We are dealing with a node, we may need to expand all dimensions
                for (int d2=0; d2<PhysicalDim; d2++) {
                    if ( range2[2*d2+1]==d_size[d2] && !d_isPeriodic[d2] && d!=PhysicalDim )
                        range2[2*d2+1]++;
                }
            } else {
                if ( d==1 && PhysicalDim==3 ) {
                    // We are dealing with a edge in 3d
                    for (int d2=0; d2<PhysicalDim; d2++) {
                        if ( range2[2*d2+1]==d_size[d2] && !d_isPeriodic[d2] && d!=PhysicalDim && s!=d2 )
                            range2[2*d2+1]++;
                    }
                } else if ( PhysicalDim-d==1 ) { 
                    // We are dealing with a face or edge in 3d or 2d respectivly
                    int d2 = s;
                    if ( range2[2*d2+1]==d_size[d2] && !d_isPeriodic[d2] && d!=PhysicalDim )
                        range2[2*d2+1]++;
                } else {
                    AMP_ERROR("Internal error");
                }
            }
            // Create the elements
            if ( PhysicalDim==3 ) {
                for (int k=range2[4]; k<range2[5]; k++) {
                    for (int j=range2[2]; j<range2[3]; j++) {
                        for (int i=range2[0]; i<range2[1]; i++) {
                            MeshElementIndex index;
                            index.type = (GeomType) d;
                            index.index[0] = i;
                            index.index[1] = j;
                            index.index[2] = k;
                            index.side = s;
                            d_elements[d][0]->push_back( index );
                        }
                    }
                }
            }
        }
    }
    if ( PhysicalDim==1 ) {
        if ( d_isPeriodic[0] )
            AMP_ASSERT( (int)d_elements[0][0]->size() == d_size[0] );
        else
            AMP_ASSERT( (int)d_elements[0][0]->size() == (d_size[0]+1) );
        AMP_ASSERT( d_comm.sumReduce((int)d_elements[1][0]->size()) == d_size[0] );
    } else if ( PhysicalDim==2 ) {
        size_t N_faces_global = d_size[0]*d_size[1];
        size_t N_edges_global = 2*d_size[0]*d_size[1];
        if ( !d_isPeriodic[0] )
            N_edges_global += d_size[1];
        if ( !d_isPeriodic[1] )
            N_edges_global += d_size[0];
        size_t N_nodes_global = 1;
        for (int i=0; i<2; i++) {
            if ( d_isPeriodic[i] )
                N_nodes_global *= d_size[i];
            else
                N_nodes_global *= d_size[i]+1;
        }
        AMP_ASSERT( d_comm.sumReduce(d_elements[0][0]->size()) == N_nodes_global );
        AMP_ASSERT( d_comm.sumReduce(d_elements[1][0]->size()) == N_edges_global );
        AMP_ASSERT( d_comm.sumReduce(d_elements[2][0]->size()) == N_faces_global );
    } else if ( PhysicalDim==3 ) {
        size_t N_elements_global = d_size[0]*d_size[1]*d_size[2];
        size_t N_faces_global = 3*d_size[0]*d_size[1]*d_size[2];
        size_t N_edges_global = 3*d_size[0]*d_size[1]*d_size[2];
        if ( !d_isPeriodic[0] ) {
            N_faces_global += d_size[1]*d_size[2];
            N_edges_global += 2*d_size[1]*d_size[2];
        }
        if ( !d_isPeriodic[1] ) {
            N_faces_global += d_size[0]*d_size[2];
            N_edges_global += 2*d_size[0]*d_size[2];
        }
        if ( !d_isPeriodic[2] ) {
            N_faces_global += d_size[0]*d_size[1];
            N_edges_global += 2*d_size[0]*d_size[1];
        }
        if ( !d_isPeriodic[0] && !d_isPeriodic[1] )
            N_edges_global += d_size[2];
        if ( !d_isPeriodic[0] && !d_isPeriodic[2] )
            N_edges_global += d_size[1];
        if ( !d_isPeriodic[1] && !d_isPeriodic[2] )
            N_edges_global += d_size[0];
        size_t N_nodes_global = 1;
        for (int i=0; i<3; i++) {
            if ( d_isPeriodic[i] )
                N_nodes_global *= d_size[i];
            else
                N_nodes_global *= d_size[i]+1;
        }
        AMP_ASSERT( d_comm.sumReduce(d_elements[0][0]->size()) == N_nodes_global );
        AMP_ASSERT( d_comm.sumReduce(d_elements[1][0]->size()) == N_edges_global );
        AMP_ASSERT( d_comm.sumReduce(d_elements[2][0]->size()) == N_faces_global );
        AMP_ASSERT( d_comm.sumReduce(d_elements[3][0]->size()) == N_elements_global );
    } else {
        AMP_ERROR("Not programmed for this dimension yet");
    }
    PROFILE_STOP("create_owned_elements");
    // Sort the elements for easy searching
    for (int d=0; d<PhysicalDim; d++)
        AMP::Utilities::quicksort( *d_elements[d][0] );
    // Create the ghost elements of type GeomType == PhysicalDim
    PROFILE_START("create_ghost_elements: 1");
    for (int gcw=1; gcw<=d_max_gcw; gcw++) {
        d_elements[PhysicalDim][gcw] = ElementIndexList( new std::vector<MeshElementIndex>() );
        if ( PhysicalDim==3 ) {
            for (int k=range[4]-gcw; k<range[5]+gcw; k++) {
                for (int j=range[2]-gcw; j<range[3]+gcw; j++) {
                    for (int i=range[0]-gcw; i<range[1]+gcw; i++) {
                        if ( ( i>range[0]-gcw && i<range[1]+gcw-1 ) &&
                             ( j>range[2]-gcw && j<range[3]+gcw-1 ) &&
                             ( k>range[4]-gcw && k<range[5]+gcw-1 ) ) {
                            // The element was already included by another ghost (or owned) cell
                            continue;
                        }
                        if ( ( (i<0||i>=d_size[0]) && !d_isPeriodic[0] ) ||
                             ( (j<0||j>=d_size[1]) && !d_isPeriodic[1] ) ||  
                             ( (k<0||k>=d_size[2]) && !d_isPeriodic[2] ) ) {
                            // The element is outside the domain
                            continue;
                        }
                        // Create the index (adjusting for periodic boundaries)
                        MeshElementIndex index( PhysicalDim, 0, i, j, k );
                        if ( i<0 ) { index.index[0] += d_size[0]; }
                        if ( j<0 ) { index.index[1] += d_size[1]; }
                        if ( k<0 ) { index.index[2] += d_size[2]; }
                        if ( i>=d_size[0] ) { index.index[0] -= d_size[0]; }
                        if ( j>=d_size[1] ) { index.index[1] -= d_size[1]; }
                        if ( k>=d_size[2] ) { index.index[2] -= d_size[2]; }
                        // Check if the element is already in one of the lists
                        bool found = false;
                        for (int k=0; k<gcw; k++) {
                            if ( d_elements[PhysicalDim][k]->size()==0 )
                                continue;
                            size_t m = AMP::Utilities::findfirst( *d_elements[PhysicalDim][k], index );
                            if ( m==d_elements[PhysicalDim][k]->size() ) { m--; }
                            if ( d_elements[PhysicalDim][k]->operator[](m) == index )
                                found = true;
                        }
                        if ( !found )
                            d_elements[PhysicalDim][gcw]->push_back( index );
                    }
                }
            }
        } else { 
            AMP_ERROR("Not programmed for this dimension yet");
        }
        // Sort the elements for easy searching and remove any duplicates
        AMP::Utilities::unique( *d_elements[PhysicalDim][gcw] );
    }
    PROFILE_STOP("create_ghost_elements: 1");
    // Create the remaining ghost elements
    PROFILE_START("create_ghost_elements: 2");
    for (int gcw=1; gcw<=d_max_gcw; gcw++) {
        // Get an iterator over all elements of the given gcw
        AMP::Mesh::MeshIterator iterator = this->getIterator( (GeomType) PhysicalDim, gcw );
        // Create the vectors that will store the ghost elements
        size_t N_elem_est = iterator.size() - numLocalElements((GeomType)PhysicalDim) - numGhostElements((GeomType)PhysicalDim,gcw-1);
        for (int d=0; d<PhysicalDim; d++) {
            d_elements[d][gcw] = ElementIndexList( new std::vector<MeshElementIndex>() );
            if ( d==0 )
                d_elements[d][gcw]->reserve( N_elem_est );
            else if ( d==1 )
                d_elements[d][gcw]->reserve( 6*N_elem_est );
            else if ( d==1 )
                d_elements[d][gcw]->reserve( 3*N_elem_est );
        }
        // Loop through the elements creating the ghosts
        for (size_t i=0; i<iterator.size(); i++) {
            // Skip any cells we know do not have ghost sub-elements
            bool interior_element = true;
            structuredMeshElement* elem = dynamic_cast<structuredMeshElement*>( iterator->getRawElement() );
            AMP_ASSERT(elem!=NULL); 
            MeshElementIndex index = elem->d_index;
            for (int d=0; d<PhysicalDim; d++) {
                // If we are outside or on the +boundary, we are not an interior element 
                if ( index.index[d]<range[2*d+0] || index.index[d]>=range[2*d+1]-1 )
                    interior_element = false;
            }
            if ( interior_element ) {
                ++iterator;
                continue;
            }
            // Check all of the sub elements to see if they may be ghosts
            for (int d=0; d<PhysicalDim; d++) {
                // Get the elements of the given type that compose the current element
                std::vector<MeshElement> elements = iterator->getElements( (GeomType) d );
                // Loop through the current elements
                for (size_t j=0; j<elements.size(); j++) {
                    elem = dynamic_cast<structuredMeshElement*>( elements[j].getRawElement() );
                    AMP_ASSERT(elem!=NULL); 
                    index = elem->d_index;
                    // Check if the current element exists in the list of elements so far
                    bool found = false;
                    for (int k=0; k<gcw; k++) {
                        if ( d_elements[d][k]->size()==0 )
                            continue;
                        size_t m = AMP::Utilities::findfirst( *d_elements[d][k], index );
                        if ( m==d_elements[d][k]->size() ) { m--; }
                        if ( d_elements[d][k]->operator[](m) == index )
                            found = true;
                    }
                    if ( !found )
                        d_elements[d][gcw]->push_back( index );
                }
            }
            ++iterator;
        }
        // Sort the elements and remove duplicates
        for (int d=0; d<PhysicalDim; d++)
            AMP::Utilities::unique( *d_elements[d][gcw] );
    }
    PROFILE_STOP("create_ghost_elements: 2");
    // Compute the number of local, global and ghost elements
    for (int d=0; d<=PhysicalDim; d++) 
        N_global[d] = d_comm.sumReduce(d_elements[d][0]->size());
    // Create the nodes
    MeshIterator nodeIterator = getIterator(Vertex,d_max_gcw);
    for (int d=0; d<PhysicalDim; d++)
        d_coord[d] = std::vector<double>(nodeIterator.size(),0.0);
    d_index = std::vector<MeshElementIndex>(nodeIterator.size());
    for (size_t i=0; i<nodeIterator.size(); i++) {
        MeshElement* elem_ptr = nodeIterator->getRawElement();
        structuredMeshElement *element = dynamic_cast<structuredMeshElement *>( elem_ptr );
        AMP_ASSERT(element!=NULL);
        MeshElementIndex index = element->d_index;
        AMP_ASSERT(index.type==0);
        d_index[i] = index;
        ++nodeIterator;
    }
    AMP::Utilities::quicksort(d_index);
    double range2[6] = { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };
    fillCartesianNodes( PhysicalDim, &d_size[0], range2, d_index, d_coord );
    // Create the list of elements on the surface
    PROFILE_START("create_surface_elements");
    for (int d=0; d<=PhysicalDim; d++) {
        d_surface_list[d] = std::vector<ElementIndexList>(d_max_gcw+1);
        for (int gcw=0; gcw<=d_max_gcw; gcw++) {
            d_surface_list[d][gcw] = boost::shared_ptr<std::vector<MeshElementIndex> >(
                new std::vector<MeshElementIndex>() );
            d_surface_list[d][gcw]->reserve( d_elements[d][gcw]->size() );
            for (size_t k=0; k<d_elements[d][gcw]->size(); k++) {
                BoxMesh::MeshElementIndex index = d_elements[d][gcw]->operator[](k);
                structuredMeshElement elem = structuredMeshElement( index, this );
                if ( elem.isOnSurface() )
                    d_surface_list[d][gcw]->push_back( index );
            }
        }
    }
    PROFILE_STOP("create_surface_elements");
    // Create the initial boundary info 
    PROFILE_START("create_boundary_elements");
    for (int side=0; side<2*PhysicalDim; side++) {
        if ( d_isPeriodic[side/2] )
            continue;
        d_ids.push_back( side );
        int R1 = 0;
        int R2 = 0;
        if ( side%2==1 ) {
            R1 = d_size[side/2]-1;
            R2 = R1+1;
        }
        for (int d=0; d<=PhysicalDim; d++) {
            std::pair<int,GeomType> id(side,(GeomType)d);
            std::vector<ElementIndexList> ghost_list(d_max_gcw+1);
            for (int gcw=0; gcw<=d_max_gcw; gcw++) {
                std::vector<MeshElementIndex> list;
                list.reserve( d_elements[d][gcw]->size() );
                for (size_t i=0; i<d_elements[d][gcw]->size(); i++) {
                    MeshElementIndex index = d_elements[d][gcw]->operator[](i);
                    if ( d==PhysicalDim ) {
                        if ( index.index[side/2]==R1 )
                            list.push_back( index );
                    } else if ( index.type==Vertex ) {
                        if ( index.index[side/2]==R2 )
                            list.push_back( index );
                    } else if ( index.type==Edge ) {
                        if ( PhysicalDim==2 ) {
                            if ( index.side==side/2 && index.index[side/2]==R2 )
                                list.push_back( index );
                        } else {
                            if ( index.side!=side/2 && index.index[side/2]==R2 )
                                list.push_back( index );
                        }
                    } else if ( index.type==Face ) {
                        if ( index.side==side/2 && index.index[side/2]==R2 )
                            list.push_back( index );
                    } else { 
                        AMP_ERROR("Unknown element type");
                    }
                }
                // Sort for searching later
                AMP::Utilities::quicksort( list );
                // Copy the list to reduce final memory usage
                ghost_list[gcw] = ElementIndexList( new std::vector<MeshElementIndex>() );
                ghost_list[gcw]->operator=( list );
            }
            std::pair<std::pair<int,GeomType>,std::vector<ElementIndexList> >  tmp( id, ghost_list );
            d_id_list.insert( tmp );
        }
    }
    PROFILE_STOP("create_boundary_elements");
    PROFILE_STOP("initialize");
}


/****************************************************************
* De-constructor                                                *
****************************************************************/
BoxMesh::~BoxMesh()
{
}


/****************************************************************
* Estimate the mesh size                                        *
****************************************************************/
size_t BoxMesh::estimateMeshSize( const MeshParameters::shared_ptr &params )
{
    // Check for valid inputs
    AMP_INSIST(params.get(),"Params must not be null");
    boost::shared_ptr<AMP::Database> db = params->getDatabase( );
    AMP_INSIST(db.get(),"Database must exist");
    size_t N_elements = 0;
    if ( db->keyExists("NumberOfElements") ) {
        // User specified the number of elements, this should override everything
        N_elements = (size_t) db->getInteger("NumberOfElements");
    } else {
        // Use the generator to determine the number of input elements
        AMP_INSIST(db->keyExists("Generator"),"Field 'Generator' must exist in database'");
        AMP_INSIST(db->keyExists("dim"),"Field 'dim' must exist in database'");
        AMP_INSIST(db->keyExists("Size"),"Field 'Size' must exist in database'");
        std::string generator = db->getString("Generator");
        int dim = db->getInteger("dim");
        std::vector<int> size = db->getIntegerArray("Size");
        for (size_t d=0; d<size.size(); d++)
            AMP_INSIST(size[d]>0,"All values of size must be > 0");
        if ( generator.compare("cube")==0 ) {
            AMP_INSIST((int)size.size()==dim,"Size of field 'Size' must match dimfor cube");
            N_elements = 1;
            for (int d=0; d<dim; d++)
                N_elements *= size[d];
        } else if ( generator.compare("tube")==0 ) {
            AMP_INSIST(dim==3,"tube requires a 3d mesh");
            AMP_INSIST((int)size.size()==dim,"Size of field 'Size' must be 1x3 for tube");
            N_elements = size[0]*size[1]*size[2];
        } else if ( generator.compare("cylinder")==0 ) {
            AMP_INSIST(dim==3,"cylinder requires a 3d mesh");
            AMP_INSIST((int)size.size()==2,"Size of field 'Size' must be 1x2 for cylinder");
            N_elements = (2*size[0])*(2*size[0])*size[1];
        } else {
            AMP_ERROR("Unkown generator");
        }
    }
    // Adjust the number of elements by a weight if desired
    if ( db->keyExists("Weight") ) {
        double weight = db->getDouble("Weight");
        N_elements = (size_t) ceil(weight*((double)N_elements));
    }
    return N_elements;
}


/****************************************************************
* Function to return the element given an ID                    *
****************************************************************/
MeshElement BoxMesh::getElement ( const MeshElementID &elem_id ) const
{
    std::vector<int> range = getLocalBlock( elem_id.owner_rank() );
    AMP_ASSERT(PhysicalDim<=3);
    // Increase the index range for the boxes on the boundary for all elements except the current dimension
    if ( elem_id.type() != PhysicalDim ) {
        for (int d=0; d<PhysicalDim; d++) {
            if ( range[2*d+1]==d_size[d] && !d_isPeriodic[d] )
                range[2*d+1]++;
        }
    }
    // Get the 3-index from the local id
    size_t myBoxSize[3]={1,1,1};
    for (int d=0; d<PhysicalDim; d++)
        myBoxSize[d] = range[2*d+1] - range[2*d+0];
    MeshElementIndex index;
    index.type = elem_id.type();
    size_t local_id = elem_id.local_id();
    index.index[0] = range[0] + (int) local_id%myBoxSize[0];
    index.index[1] = range[2] + (int) (local_id/myBoxSize[0])%myBoxSize[1];
    index.index[2] = range[4] + (int) (local_id/(myBoxSize[0]*myBoxSize[1]))%myBoxSize[2];
    index.side = (unsigned char) (local_id/(myBoxSize[0]*myBoxSize[1]*myBoxSize[2]));
    // Create the element
    structuredMeshElement elem( index, this );
    AMP_ASSERT(elem.globalID()==elem_id);
    return elem;
}


/****************************************************************
* Functions to return the number of elements                    *
****************************************************************/
size_t BoxMesh::numLocalElements( const GeomType type ) const
{
    return d_elements[(int)type][0]->size();
}
size_t BoxMesh::numGlobalElements( const GeomType type ) const
{
    return N_global[(int)type];
}
size_t BoxMesh::numGhostElements( const GeomType type, int gcw ) const
{
    size_t N_ghost = 0;
    for (int i=1; i<=gcw; i++)
        N_ghost += d_elements[(int)type][gcw]->size();
    return N_ghost;
}


/****************************************************************
* Function to get an iterator                                   *
****************************************************************/
MeshIterator BoxMesh::getIterator( const GeomType type, const int gcw ) const
{
    AMP_ASSERT(type<=3);
    AMP_ASSERT(gcw<(int)d_elements[type].size());
    // Construct a list of iterators over the elements of interest
    std::vector<boost::shared_ptr<MeshIterator> > iterator_list;
    iterator_list.reserve(gcw+1);
    for (int i=0; i<=gcw; i++) {
        iterator_list.push_back( boost::shared_ptr<MeshIterator>(
            new structuredMeshIterator( d_elements[type][i], this, 0 ) ) );
    }
    // Create the iterator
    MultiIterator iterator(iterator_list, 0 );
    return iterator;
}


/****************************************************************
* Function to get an iterator over the surface                  *
****************************************************************/
MeshIterator BoxMesh::getSurfaceIterator( const GeomType type, const int gcw ) const
{
    AMP_ASSERT(type<=3);
    AMP_ASSERT(gcw<(int)d_surface_list[type].size());
    // Construct a list of iterators over the elements of interest
    std::vector<boost::shared_ptr<MeshIterator> > iterator_list;
    iterator_list.reserve(gcw+1);
    for (int i=0; i<=gcw; i++) {
        iterator_list.push_back( boost::shared_ptr<MeshIterator>(
            new structuredMeshIterator( d_surface_list[type][i], this, 0 ) ) );
    }
    // Create the iterator
    return MultiIterator(iterator_list, 0 );
}


/****************************************************************
* Functions to get the boundaries                               *
****************************************************************/
std::vector<int> BoxMesh::getBoundaryIDs ( ) const
{
    return d_ids;
}
MeshIterator BoxMesh::getBoundaryIDIterator ( const GeomType type, const int id, const int gcw) const
{
    std::map<std::pair<int,GeomType>,std::vector<ElementIndexList> >::const_iterator it = d_id_list.find( std::pair<int,GeomType>(id,type) );
    if ( it==d_id_list.end() )
        return MeshIterator();
    AMP_ASSERT(gcw<(int)it->second.size());
    // Construct a list of iterators over the elements of interest
    std::vector<boost::shared_ptr<MeshIterator> > iterator_list;
    iterator_list.reserve(gcw+1);
    for (int i=0; i<=gcw; i++) {
        iterator_list.push_back( boost::shared_ptr<MeshIterator>(
            new structuredMeshIterator( it->second[i], this, 0 ) ) );
    }
    // Create the iterator
    return MultiIterator(iterator_list, 0 );
}
std::vector<int> BoxMesh::getBlockIDs ( ) const
{
    return std::vector<int>(1,0);
}
MeshIterator BoxMesh::getBlockIDIterator ( const GeomType type, const int id, const int gcw ) const
{
    if ( id==0 ) 
        return getIterator( type, gcw );
    return MeshIterator();
}


/****************************************************************
* Functions to displace the mesh                                *
****************************************************************/
void BoxMesh::displaceMesh( std::vector<double> x )
{
    AMP_ASSERT(x.size()==PhysicalDim);
    for (int i=0; i<PhysicalDim; i++) {
        for (size_t j=0; j<d_coord[i].size(); j++)
            d_coord[i][j] += x[i];
        d_box[2*i+0] += x[i];
        d_box[2*i+1] += x[i];
        d_box_local[2*i+0] += x[i];
        d_box_local[2*i+1] += x[i];
    }
}
#ifdef USE_AMP_VECTORS
void BoxMesh::displaceMesh( const AMP::LinearAlgebra::Vector::const_shared_ptr x )
{
    // Create the position vector with the necessary ghost nodes
    AMP::Discretization::DOFManager::shared_ptr DOFs = 
    AMP::Discretization::simpleDOFManager::create( 
        shared_from_this(), getIterator(AMP::Mesh::Vertex,d_max_gcw), getIterator(AMP::Mesh::Vertex,0), PhysicalDim );
    AMP::LinearAlgebra::Variable::shared_ptr nodalVariable( new AMP::LinearAlgebra::Variable( "tmp_pos" ) );
    AMP::LinearAlgebra::Vector::shared_ptr displacement = AMP::LinearAlgebra::createVector( DOFs, nodalVariable, false );
    std::vector<size_t> dofs1(PhysicalDim);
    std::vector<size_t> dofs2(PhysicalDim);
    AMP::Mesh::MeshIterator cur = getIterator(AMP::Mesh::Vertex,0);
    AMP::Mesh::MeshIterator end = cur.end();
    AMP::Discretization::DOFManager::shared_ptr DOFx = x->getDOFManager();
    std::vector<double> data(PhysicalDim);
    while ( cur != end ) {
        AMP::Mesh::MeshElementID id = cur->globalID();
        DOFx->getDOFs( id, dofs1 );
        DOFs->getDOFs( id, dofs2 );
        x->getValuesByGlobalID( PhysicalDim, &dofs1[0], &data[0] );
        displacement->setValuesByGlobalID( PhysicalDim, &dofs2[0], &data[0] );
        ++cur;
    }
    displacement->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    // Move all nodes (including the ghost nodes)
    std::vector<size_t> dofs(PhysicalDim);
    std::vector<double> disp(PhysicalDim);
    for (size_t i=0; i<d_coord[0].size(); i++) {
        MeshElementID id = structuredMeshElement( d_index[i], this ).globalID();
        DOFs->getDOFs( id, dofs );
        AMP_ASSERT(dofs.size()==PhysicalDim);
        displacement->getValuesByGlobalID( (int)PhysicalDim, &dofs[0], &disp[0] );
        for (int j=0; j<PhysicalDim; j++)
            d_coord[j][i] += disp[j];
    }
    // Compute the new bounding box of the mesh
    d_box_local = std::vector<double>(2*PhysicalDim);
    for (int d=0; d<PhysicalDim; d++) {
        d_box_local[2*d+0] = 1e100;
        d_box_local[2*d+1] = -1e100;
        for (size_t i=0; i<d_coord[d].size(); i++) {
            if ( d_coord[d][i]<d_box_local[2*d+0])
                d_box_local[2*d+0] = d_coord[d][i];
            if ( d_coord[d][i]>d_box_local[2*d+1])
                d_box_local[2*d+1] = d_coord[d][i];
        }
    }
    d_box = std::vector<double>(PhysicalDim*2);
    for (int i=0; i<PhysicalDim; i++) {
        d_box[2*i+0] = d_comm.minReduce( d_box_local[2*i+0] );
        d_box[2*i+1] = d_comm.maxReduce( d_box_local[2*i+1] );
    } 
}
#endif


/****************************************************************
* Helper function to return the indices of the local block      *
* owned by the given processor                                  *
****************************************************************/
std::vector<int> BoxMesh::getLocalBlock(unsigned int rank) const
{
    size_t num_blocks = 1;
    for (int d=0; d<PhysicalDim; d++)
        num_blocks *= d_numBlocks[d];
    AMP_ASSERT((int)rank<(int)num_blocks);
    std::vector<int> range(2*PhysicalDim);
    int tmp = 1;
    for (int d=0; d<PhysicalDim; d++) {
        size_t i = (size_t) ((rank/tmp)%d_numBlocks[d]);
        tmp *= d_numBlocks[d];
        size_t size = (size_t) d_size[d];
        size_t N_blocks = (size_t) d_numBlocks[d];
        range[2*d+0] = (int) ((i*size)/N_blocks);
        range[2*d+1] = (int) (((i+1)*size)/N_blocks);
        range[2*d+1] = std::min(range[2*d+1],d_size[d]);
    }
    return range;
}


/****************************************************************
* Helper function to return the indices and rank of the owning  *
* block for a given MeshElementIndex                            *
****************************************************************/
void BoxMesh::getOwnerBlock(const MeshElementIndex index, unsigned int &rank, int *range ) const
{
    int myBoxIndex[3]={1,1,1};
    for (int d=0; d<PhysicalDim; d++) {
        size_t size = (size_t) d_size[d];
        size_t N_blocks = (size_t) d_numBlocks[d];
        if ( index.index[d]==d_size[d] ) {
            // The element lies on the physical bounadry
            AMP_ASSERT(index.type<PhysicalDim);
            myBoxIndex[d] = d_numBlocks[d]-1;
        } else {
            // Find the owning box
            myBoxIndex[d] = (int) ((((size_t)index.index[d]+1)*N_blocks-1)/size);
        }
        range[2*d+0] = (int) ((size*((size_t)myBoxIndex[d]))/N_blocks);
        range[2*d+1] = (int) ((size*((size_t)myBoxIndex[d]+1))/N_blocks);
        range[2*d+1] = std::min(range[2*d+1],d_size[d]);
    }
    // Increase the index range for the boxes on the boundary for all elements except the current dimension
    if ( index.type != PhysicalDim ) {
        for (int d=0; d<PhysicalDim; d++) {
            if ( range[2*d+1]==d_size[d] && !d_isPeriodic[d] )
                range[2*d+1]++;
        }
    }
    rank = (unsigned int) ( myBoxIndex[0] + myBoxIndex[1]*d_numBlocks[0] + 
        myBoxIndex[2]*d_numBlocks[0]*d_numBlocks[1] );
}


/****************************************************************
* Helper function to fill the cartesian coordinates             *
****************************************************************/
void BoxMesh::fillCartesianNodes(int dim, const int* globalSize, const double *range, 
    const std::vector<MeshElementIndex> &index, std::vector<double> *coord)
{
    AMP_ASSERT(index.size()==coord[0].size());
    for (size_t i=0; i<index.size(); i++) {
        AMP_ASSERT(index[i].type==0);
        for (int d=0; d<dim; d++)
            coord[d][i] = range[2*d+0] + (range[2*d+1]-range[2*d+0])*((double)index[i].index[d])/((double)globalSize[d]);
    }
}


} // Mesh namespace
} // AMP namespace

