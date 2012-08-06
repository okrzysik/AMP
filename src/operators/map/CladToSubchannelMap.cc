#if 1
#include "operators/map/CladToSubchannelMap.h"
#include "discretization/DOF_Manager.h"
#include "utils/PIO.h"
#include "utils/ProfilerApp.h"


namespace AMP {
namespace Operator {

static double interp_linear(const std::vector<double>&, const std::vector<double>&, double );


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

    // Clone the communicator to protect the communication (we need a large number of unique tags)
    d_MapComm = d_MapComm.dup();
    d_currRequests = std::vector<MPI_Request>();

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

    // Create the send/recv buffers
    d_sendBuffer = std::vector<std::vector<std::pair<double,int> > >(N_subchannels);
    if ( d_mesh1.get() != NULL ) {
        for (size_t i=0; i<N_subchannels; i++) {
            if ( d_elem[i].size() > 0 )
                d_sendBuffer[i] = std::vector<std::pair<double,int> >(d_z.size());
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
    d_ownSubChannel = std::vector<bool>(Nx*Ny,false);
    if ( d_mesh2.get() != NULL ) {
        AMP::Mesh::MeshIterator it = d_iterator2.begin();
        AMP_ASSERT(it.size()>0);
        for (size_t k=0; k<it.size(); k++) {
            std::vector<double> center = it->centroid();
            size_t i = Utilities::findfirst(d_x,center[0]);
            size_t j = Utilities::findfirst(d_y,center[1]);
            AMP_ASSERT( i>0 && i<d_x.size() && j>0 && j<d_y.size() );
            d_ownSubChannel[(i-1)+(j-1)*Nx] = true;
            ++it;
        }
    }
    d_subchannelRanks = std::vector<std::vector<int> >(Nx*Ny);
    std::vector<char> tmp(d_MapComm.getSize());
    for (size_t i=0; i<Nx*Ny; i++) {
        d_MapComm.allGather((char)(d_ownSubChannel[i]?1:0),&tmp[0]);
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
void CladToSubchannelMap::applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &,
    const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &,
    const double, const double)
{
    // Fill the send buffer
    if ( d_mesh1.get() != NULL ) {
        AMP::Discretization::DOFManager::shared_ptr  DOF = u->getDOFManager( );
        std::vector<size_t> dofs;
        for (size_t i=0; i<N_subchannels; i++) {
            if ( d_subchannelSend.empty() )
                continue;
            for (size_t j=0; j<d_z.size(); j++)
                d_sendBuffer[i][j] = std::pair<double,int>(0.0,0);
            for (size_t j=0; j<d_elem[i].size(); j++) {
                DOF->getDOFs(d_elem[i][j],dofs);
                AMP_ASSERT(dofs.size()==1);
                std::vector<double> coord = d_mesh1->getElement(d_elem[i][j]).centroid();
                size_t index = AMP::Utilities::findfirst(d_z,coord[2]-1e-10);
                AMP_ASSERT(index<d_z.size());
                AMP_ASSERT(Utilities::approx_equal(d_z[index],coord[2]));
                double val = u->getLocalValueByGlobalID(dofs[0]);
                d_sendBuffer[i][index].first += val;
                d_sendBuffer[i][index].second++;
            }
        }
    }
    // Send the data
    if ( d_mesh1.get() != NULL ) {
        for (size_t i=0; i<N_subchannels; i++) {
            if ( d_elem[i].empty() )
                continue;
            int tag = (int) i;  // We have an independent comm
            for (size_t j=0; j<d_subchannelRanks[i].size(); j++) {
                int rank = d_subchannelRanks[i][j];
                d_currRequests.push_back( d_MapComm.Isend( &d_sendBuffer[i], (int)d_z.size(),rank, tag ) );
            }
        }
    }    
}


/************************************************************************
*  Finish the communication                                             *
************************************************************************/
void CladToSubchannelMap::applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &,
    const AMP::LinearAlgebra::Vector::shared_ptr &, AMP::LinearAlgebra::Vector::shared_ptr &,
    const double, const double)
{
    // Recieve the data
    int length = (int) d_z.size();
    std::vector<std::vector<double> > scalarZmaps(N_subchannels);
    std::vector<std::pair<double,int> > data1, data2;
    for (size_t i=0; i<N_subchannels; i++) {
        if ( d_ownSubChannel[i] ) {
            int tag = (int) i;  // We have an independent comm
            for (size_t k=0; k<d_z.size(); k++) {
                data1[k].first = 0.0;
                data1[k].second = 0;
            }
            for (size_t j=0; j<d_subchannelSend[i].size(); j++) {
                d_MapComm.recv( &data2[0],length, d_subchannelSend[i][j], false, tag );
                for (size_t k=0; k<d_z.size(); k++) {
                    data1[k].first += data2[k].first;
                    data1[k].second += data2[k].second;
                }
            }
            for (size_t k=0; k<d_z.size(); k++) {
                AMP_ASSERT(data1[k].second>0);
                scalarZmaps[i][k] = data1[k].first/data1[k].second;
            }
        }
    }
    // Fill the output vector
    AMP::Discretization::DOFManager::shared_ptr  DOF = d_OutputVector->getDOFManager( );
    AMP::Mesh::MeshIterator it = d_iterator2.begin();
    std::vector<size_t> dofs;
    for (size_t i=0; i<it.size(); i++) {
        DOF->getDOFs(it->globalID(),dofs);
        AMP_ASSERT(dofs.size()==1);
        std::vector<double> pos = it->centroid();
        double val = interp_linear(d_z,scalarZmaps[i],pos[2]);
        d_OutputVector->setLocalValueByGlobalID(dofs[0],val);
        ++it;
    }
    // Wait for all communication to finish
    if ( d_currRequests.size() > 0 )
        AMP::AMP_MPI::waitAll( (int)d_currRequests.size(), &d_currRequests[0] );
    d_currRequests.resize(0);
    // Call makeConsistent
    d_OutputVector->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
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


/************************************************************************
*  1D interpolation                                                     *
************************************************************************/
double interp_linear(const std::vector<double>& x, const std::vector<double>& f, double pos )
{
    AMP_ASSERT(x.size()>=2&&x.size()==f.size());
    size_t i = AMP::Utilities::findfirst(x,pos);
    if ( i==0 ) { i = 1; }
    return f[i-1] + (pos-x[i-1])/(x[i]-x[i-1])*f[i];
}



} // Operator namespace
} // AMP namespace

#endif
