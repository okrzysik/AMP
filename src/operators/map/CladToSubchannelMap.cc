#if 1
#include "operators/map/CladToSubchannelMap.h"
#include "discretization/DOF_Manager.h"
#include "utils/PIO.h"
#include "utils/ProfilerApp.h"


namespace AMP {
namespace Operator {


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
CladToSubchannelMap::CladToSubchannelMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &p )
    : AsyncMapOperator ( p )
{
    boost::shared_ptr<CladToSubchannelMapParameters>  params = boost::dynamic_pointer_cast<CladToSubchannelMapParameters> ( p );
    AMP_ASSERT( params );

    int DofsPerObj = params->d_db->getInteger ( "DOFsPerObject" );
    AMP_INSIST(DofsPerObj==1,"CladToSubchannelMap is currently only designed for 1 DOF per node");

    // Get the iterators
    if ( d_mesh1.get() != NULL )
        d_iterator1 = d_mesh1->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID1, 0 );
    if ( d_mesh2.get() != NULL )
        d_iterator2 = getSubchannelIterator();
    
    // Get the x-y-z grid for the subchannel mesh
    fillSubchannelGrid();
    
    // For each subchannel, get the list of local MeshElement in that channel
    size_t Nx = d_x.size()-1;
    d_elem = std::vector<std::vector<AMP::Mesh::MeshElementID> >(N_subchannels);
    if ( d_mesh1.get() != NULL ) {
        AMP::Mesh::MeshIterator it = d_iterator1.begin();
        for (size_t k=0; k<it.size(); k++) {
            std::vector<double> center = it->centroid();
            size_t i = Utilities::findfirst(d_x,center[0]);
            size_t j = Utilities::findfirst(d_y,center[1]);
            if ( i>0 && i<d_x.size() && j>0 && j<d_y.size() )
                d_elem[(i-1)+(j-1)*Nx].push_back( it->globalID() );
            ++it;
        }
    }
    
    // Get the list of processors that we will recieve from for each subchannel
    d_subchannelSend = std::vector<std::vector<int> >(N_subchannels);
    std::vector<char> tmp(d_MapComm.getSize());
    for (size_t i=0; i<N_subchannels; i++) {
        d_MapComm.allGather((char)(d_elem[i].size()>0?1:0),&tmp[0]);
        for (size_t j=0; j<tmp.size(); j++) {
            if ( tmp[j]==1 )
                d_subchannelSend[i].push_back(j);
        }
    }
    
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
CladToSubchannelMap::~CladToSubchannelMap ()
{
}


/************************************************************************
*  Check if the map type is "CladToSubchannelMap"                       *
************************************************************************/
bool CladToSubchannelMap::validMapType ( const std::string &t )
{
    if ( t == "CladToSubchannelMap" )
        return true;
    return false;
}


/************************************************************************
*  Function to fill the grid of the subchannel for all processors       *
************************************************************************/
void CladToSubchannelMap::fillSubchannelGrid()
{
    // Create the grid for all processors
    std::set<double> x, y, z;
    if ( d_mesh2.get() != NULL ) {
        AMP::Mesh::MeshIterator it = d_mesh2->getIterator( AMP::Mesh::Vertex, 0 );
        for (size_t i=0; i<it.size(); i++) {
            std::vector<double> coord = it->coord();
            AMP_ASSERT(coord.size()==3);
            x.insert( coord[0] );
            y.insert( coord[1] );
            z.insert( coord[2] );
            ++it;
        }
    }
    d_MapComm.setGather(x);
    d_MapComm.setGather(y);
    d_MapComm.setGather(z);
    double last = 1e300;
    for (std::set<double>::iterator it=x.begin(); it!=x.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            x.erase(it);
        else
            last = *it;
    }
    for (std::set<double>::iterator it=y.begin(); it!=y.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            y.erase(it);
        else
            last = *it;
    }
    for (std::set<double>::iterator it=z.begin(); it!=z.end(); ++it) {
        if ( Utilities::approx_equal(last,*it,1e-12) )
            z.erase(it);
        else
            last = *it;
    }
    d_x = std::vector<double>(x.begin(),x.end());
    d_y = std::vector<double>(y.begin(),y.end());
    d_z = std::vector<double>(z.begin(),z.end());
    size_t Nx = d_x.size()-1;
    size_t Ny = d_y.size()-1;
    size_t Nz = d_z.size()-1;
    if ( d_mesh2.get() != NULL ) 
        AMP_ASSERT(Nx*Ny*Nz==d_mesh2->numGlobalElements(AMP::Mesh::Volume));
    N_subchannels = Nx*Ny;
    // Get a list of processors that need each x-y point
    std::vector<bool> own_xy(Nx*Ny,false);
    if ( d_mesh2.get() != NULL ) {
        AMP::Mesh::MeshIterator it = d_iterator2.begin();
        AMP_ASSERT(it.size()>0);
        for (size_t k=0; k<it.size(); k++) {
            std::vector<double> center = it->centroid();
            size_t i = Utilities::findfirst(d_x,center[0]);
            size_t j = Utilities::findfirst(d_y,center[1]);
            AMP_ASSERT( i>0 && i<d_x.size() && j>0 && j<d_y.size() );
            own_xy[(i-1)+(j-1)*Nx] = true;
            ++it;
        }
    }
    d_subchannelRanks = std::vector<std::vector<int> >(Nx*Ny);
    std::vector<char> tmp(d_MapComm.getSize());
    for (size_t i=0; i<Nx*Ny; i++) {
        d_MapComm.allGather((char)(own_xy[i]?1:0),&tmp[0]);
        for (size_t j=0; j<tmp.size(); j++) {
            if ( tmp[j]==1 )
                d_subchannelRanks[i].push_back(j);
        }
    }
}


/************************************************************************
*  Get an iterator over the desired faces in the subchannel mesh        *
************************************************************************/
AMP::Mesh::MeshIterator CladToSubchannelMap::getSubchannelIterator()
{
    std::multimap<double,AMP::Mesh::MeshElement> xyFace;
    AMP::Mesh::MeshIterator iterator = d_mesh2->getIterator( AMP::Mesh::Face, 0 );
    for(size_t i=0; i<iterator.size(); ++i ) {
        std::vector<AMP::Mesh::MeshElement> nodes = iterator->getElements(AMP::Mesh::Vertex);
        std::vector<double> center = iterator->centroid();
        bool is_valid = true;
        for (size_t j=0; j<nodes.size(); ++j) {
            std::vector<double> coord = nodes[j].coord();
            if ( !AMP::Utilities::approx_equal(coord[2],center[2], 1e-6) )
                is_valid = false;
        }
        if ( is_valid ) {
            xyFace.insert(std::pair<double,AMP::Mesh::MeshElement>(center[2],*iterator));
        }
        ++iterator;
    }
    boost::shared_ptr<std::vector<AMP::Mesh::MeshElement> > elements( new std::vector<AMP::Mesh::MeshElement>() );
    elements->reserve(xyFace.size());
    for (std::multimap<double,AMP::Mesh::MeshElement>::iterator it=xyFace.begin(); it!=xyFace.end(); ++it)
        elements->push_back( it->second );
    return AMP::Mesh::MultiVectorIterator( elements );
}


/************************************************************************
*  Start the communication                                              *
************************************************************************/
void CladToSubchannelMap::applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &f,
    const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
    const double a, const double b)
{
    AMP_ERROR("Not finished");
}


/************************************************************************
*  Finish the communication                                             *
************************************************************************/
void CladToSubchannelMap::applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &f,
    const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
    const double a, const double b)
{
    AMP_ERROR("Not finished");
}


/************************************************************************
*  Set the vector                                                       *
************************************************************************/
void  CladToSubchannelMap::setVector( AMP::LinearAlgebra::Vector::shared_ptr &result )
{
    if ( result.get() != NULL )
        d_OutputVector = subsetInputVector( result );
    else
        d_OutputVector = AMP::LinearAlgebra::Vector::shared_ptr();
    if ( d_mesh2.get() != NULL )
        AMP_ASSERT(d_OutputVector.get()!=NULL);
}

} // Operator namespace
} // AMP namespace

#endif
