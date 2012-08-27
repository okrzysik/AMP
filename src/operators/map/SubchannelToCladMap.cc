#if 1
#include "operators/map/SubchannelToCladMap.h"
#include "discretization/DOF_Manager.h"
#include "utils/PIO.h"
#include "utils/ProfilerApp.h"


namespace AMP {
namespace Operator {


static double interp_linear(const std::vector<double>&, const std::vector<double>&, double );


/************************************************************************
*  Default constructor                                                  *
************************************************************************/
SubchannelToCladMap::SubchannelToCladMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &p )
    : AsyncMapOperator ( p )
{
    boost::shared_ptr<SubchannelToCladMapParameters>  params = boost::dynamic_pointer_cast<SubchannelToCladMapParameters> ( p );
    AMP_ASSERT( params );

    int DofsPerObj = params->d_db->getInteger ( "DOFsPerObject" );
    AMP_INSIST(DofsPerObj==1,"SubchannelToCladMap is currently only designed for 1 DOF per node");

    // Clone the communicator to protect the communication (we need a large number of unique tags)
    d_MapComm = d_MapComm.dup();
    d_currRequests = std::vector<MPI_Request>();

    // Get the iterators
    if ( d_mesh1.get() != NULL )
        d_iterator1 = getSubchannelIterator(d_mesh1);
    if ( d_mesh2.get() != NULL )
        d_iterator2 = d_mesh2->getBoundaryIDIterator( AMP::Mesh::Vertex, params->d_BoundaryID1, 0 );
    
    // Get the x-y-z grid for the subchannel mesh
    fillSubchannelGrid(d_mesh1);
    
    // For each subchannel, get the list of local MeshElement in that channel
    d_elem = std::vector<std::vector<AMP::Mesh::MeshElementID> >(N_subchannels);
    if ( d_mesh2.get() != NULL ) {
        AMP::Mesh::MeshIterator it = d_iterator2.begin();
        for (size_t k=0; k<it.size(); k++) {
            std::vector<double> center = it->centroid();
            int index = getSubchannelIndex( center[0], center[1] );
            if ( index>=0 )
                d_elem[index].push_back( it->globalID() );
            ++it;
        }
    }
    
    // Get the list of processors that we will recieve from for each subchannel
    d_subchannelRecv = std::vector<std::vector<int> >(N_subchannels);
    std::vector<char> tmp(d_MapComm.getSize());
    for (size_t i=0; i<N_subchannels; i++) {
        d_MapComm.allGather((char)(d_elem[i].size()>0?1:0),&tmp[0]);
        for (size_t j=0; j<tmp.size(); j++) {
            if ( tmp[j]==1 )
                d_subchannelRecv[i].push_back(j);
        }
    }

    // Create the send/recv buffers
    d_sendBuffer = std::vector<std::vector<double> >(N_subchannels);
    if ( d_mesh1.get() != NULL ) {
        for (size_t i=0; i<N_subchannels; i++) {
            if ( d_ownSubChannel[i] )
                d_sendBuffer[i] = std::vector<double>(d_z.size());
        }
    }
}


/************************************************************************
*  De-constructor                                                       *
************************************************************************/
SubchannelToCladMap::~SubchannelToCladMap ()
{
}


/************************************************************************
*  Check if the map type is "SubchannelToCladMap"                       *
************************************************************************/
bool SubchannelToCladMap::validMapType ( const std::string &t )
{
    if ( t == "SubchannelToCladMap" )
        return true;
    return false;
}


/************************************************************************
*  Function to fill the grid of the subchannel for all processors       *
************************************************************************/
void SubchannelToCladMap::fillSubchannelGrid(AMP::Mesh::Mesh::shared_ptr mesh)
{
    // Create the grid for all processors
    std::set<double> x, y, z;
    if ( mesh.get() != NULL ) {
        AMP::Mesh::MeshIterator it = mesh->getIterator( AMP::Mesh::Vertex, 0 );
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
    if ( mesh.get() != NULL ) 
        AMP_ASSERT(Nx*Ny*Nz==mesh->numGlobalElements(AMP::Mesh::Volume));
    N_subchannels = Nx*Ny;
    // Get a list of processors that need each x-y point
    d_ownSubChannel = std::vector<bool>(Nx*Ny,false);
    if ( mesh.get() != NULL ) {
        AMP::Mesh::MeshIterator it = getSubchannelIterator(mesh);
        AMP_ASSERT(it.size()>0);
        for (size_t k=0; k<it.size(); k++) {
            std::vector<double> center = it->centroid();
            int index = getSubchannelIndex( center[0], center[1] );
            AMP_ASSERT(index>=0);
            d_ownSubChannel[index] = true;
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
AMP::Mesh::MeshIterator SubchannelToCladMap::getSubchannelIterator(AMP::Mesh::Mesh::shared_ptr mesh)
{
    std::multimap<double,AMP::Mesh::MeshElement> xyFace;
    AMP::Mesh::MeshIterator iterator = mesh->getIterator( AMP::Mesh::Face, 0 );
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
void SubchannelToCladMap::applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr,
    AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr,
    const double, const double)
{
    // Check if we have any data to send
    if ( d_mesh1.get() == NULL )
        return;

    // Subset the vector for the variable (we only need the local portion of the vector)
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    AMP::LinearAlgebra::VS_Comm commSelector( AMP_MPI(AMP_COMM_SELF) );
    AMP::LinearAlgebra::Vector::const_shared_ptr  commSubsetVec = u->constSelect(commSelector, var->getName());
    AMP::LinearAlgebra::Vector::const_shared_ptr  curPhysics = commSubsetVec->constSubsetVectorForVariable(var);
    AMP_ASSERT(curPhysics);

    // Zero out the send buffers
    for (size_t i=0; i<N_subchannels; i++) {
        if ( !d_ownSubChannel[i] )
            continue;
        for (size_t j=0; j<d_z.size(); j++)
            d_sendBuffer[i][j] = 0.0;
    }
    // Fill the faces I own
    AMP::Discretization::DOFManager::shared_ptr  DOF = curPhysics->getDOFManager( );
    std::vector<size_t> dofs;
    AMP::Mesh::MeshIterator it = d_iterator1.begin();
    for (size_t i=0; i<it.size(); i++) {
        AMP::Mesh::MeshElementID id = it->globalID();
        if ( !id.is_local() )
            continue;
        DOF->getDOFs(id,dofs);
        AMP_ASSERT(dofs.size()==1);
        double val = curPhysics->getLocalValueByGlobalID(dofs[0]);
        std::vector<double> pos = it->centroid();
        int index = getSubchannelIndex(pos[0],pos[1]);
        AMP_ASSERT(index>=0);
        AMP_ASSERT(d_ownSubChannel[index]);
        int index2 = Utilities::findfirst(d_z,pos[2]-1e-12);
        AMP_ASSERT(Utilities::approx_equal(pos[2],d_z[index2]));
        d_sendBuffer[index][index2] = val;
        ++it;
    }

    // Send the data
    for (size_t i=0; i<N_subchannels; i++) {
        if ( !d_ownSubChannel[i] )
            continue;
        int tag = (int) i;  // We have an independent comm
        for (size_t j=0; j<d_subchannelRecv[i].size(); j++) {
            int rank = d_subchannelRecv[i][j];
            d_currRequests.push_back( d_MapComm.Isend<double>( &d_sendBuffer[i][0], d_sendBuffer[i].size(), rank, tag ) );
        }
    }    
}


/************************************************************************
*  Finish the communication                                             *
************************************************************************/
void SubchannelToCladMap::applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr,
    AMP::LinearAlgebra::Vector::const_shared_ptr, AMP::LinearAlgebra::Vector::shared_ptr,
    const double, const double)
{
    if ( d_mesh2.get() == NULL ) {
        // We don't have an output vector to fill, wait for communication to finish and return
        if ( d_currRequests.size() > 0 )
            AMP::AMP_MPI::waitAll( (int)d_currRequests.size(), &d_currRequests[0] );
        d_currRequests.resize(0);
        return;
    }
    // Recieve the data
    std::vector<std::vector<double> >  f(N_subchannels);
    double *tmp_data = new double[d_z.size()];
    for (size_t i=0; i<N_subchannels; i++) {
        if ( d_elem[i].size() > 0 ) {
            f[i] = std::vector<double>(d_z.size(),0.0);
            int tag = (int) i;  // We have an independent comm
            for (size_t j=0; j<d_subchannelRanks[i].size(); j++) {
                int length = d_z.size();
                d_MapComm.recv<double>( tmp_data, length, d_subchannelRanks[i][j], false, tag );
                for (size_t k=0; k<d_z.size(); k++)
                    f[i][k] += tmp_data[k];     // This works since each z-point in only owned by one processor
            }
        }
    }
    delete [] tmp_data;
    // Fill the output vector
    AMP::Discretization::DOFManager::shared_ptr  DOF = d_OutputVector->getDOFManager( );
    std::vector<size_t> dofs;
    for (size_t i=0; i<N_subchannels; i++) {
        for (size_t j=0; j<d_elem[i].size(); j++) {
            DOF->getDOFs(d_elem[i][j],dofs);
            AMP_ASSERT(dofs.size()==1);
            std::vector<double> pos = d_mesh2->getElement(d_elem[i][j]).coord();
            double val = interp_linear(d_z,f[i],pos[2]);
            d_OutputVector->setLocalValueByGlobalID(dofs[0],val);
        }
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
void  SubchannelToCladMap::setVector( AMP::LinearAlgebra::Vector::shared_ptr &result )
{
    if ( result.get() != NULL )
        d_OutputVector = subsetOutputVector( result );
    else
        d_OutputVector = AMP::LinearAlgebra::Vector::shared_ptr();
    if ( d_mesh2.get() != NULL )
        AMP_ASSERT(d_OutputVector.get()!=NULL);
}


/************************************************************************
*  Misc functions                                                       *
************************************************************************/
int SubchannelToCladMap::getSubchannelIndex( double x, double y )
{
    size_t i = Utilities::findfirst(d_x,x);
    size_t j = Utilities::findfirst(d_y,y);
    if ( i>0 && i<d_x.size() && j>0 && j<d_y.size() )
        return (i-1)+(j-1)*(d_x.size()-1);
    return -1;
}
double interp_linear(const std::vector<double>& x, const std::vector<double>& f, double pos )
{
    AMP_ASSERT(x.size()>=2&&x.size()==f.size());
    size_t i = AMP::Utilities::findfirst(x,pos);
    if ( i==0 ) { i = 1; }
    return f[i-1] + (pos-x[i-1])/(x[i]-x[i-1])*(f[i]-f[i-1]);
}



} // Operator namespace
} // AMP namespace

#endif
