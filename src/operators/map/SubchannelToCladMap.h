//#ifndef included_AMP_SubchannelToCladMap
#define included_AMP_SubchannelToCladMap

#include "operators/map/ScalarZAxisMap.h"
#include "ampmesh/MeshElementVectorIterator.h"

namespace AMP {
namespace Operator {


typedef AsyncMapOperatorParameters  SubchannelToCladMapParameters;


/**
 * \class  SubchannelToCladMap
 * \brief  A class used to map temperature from the clad meshes to the subchannel mesh
 * \details  This class maps a scalar quantity (eg. Temperature) from the outer surface
 *    of the clad meshes of a fuel assembly to a subchannel mesh.
 *    mesh1 - clad meshes
 *    mesh2 - subchannel mesh
 */
class SubchannelToCladMap : public AMP::Operator::AsyncMapOperator
{
public:

    /** \brief  Returns true if MapType = "SubchannelToCladMapParameters"
     * \param[in] s  A string extracted from the MapType line in a MeshToMeshMap db
     * \return  True iff s == "SubchannelToCladMapParameters"
     */
    static bool validMapType ( const std::string &s );

    /** \brief  Typedef to identify the parameters class of this operator
     */
    typedef  SubchannelToCladMapParameters   Parameters;

    /** \brief   Standard constructor
     * \param[in] params  Input parameters
     */
    SubchannelToCladMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &params );

    //! Destructor
    ~SubchannelToCladMap();

    /** \brief   Set a frozen vector for results of the apply operation. 
     * \details  Set a frozen vector for results of the apply operation. 
     * \param result    The results vector
     */
    virtual void  setVector ( AMP::LinearAlgebra::Vector::shared_ptr &result );

    /** \brief   Start a communicative apply operation. 
     * \details  Start a communicative apply operation. 
     */
    virtual void applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a = -1.0, const double b = 1.0);

    /** \brief   Finish a communicative apply operation. 
     * \details  Finish a communicative apply operation. 
     */
    virtual void applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a = -1.0, const double b = 1.0);

protected:


private:

    // The grid of the subchannel mesh
    void fillSubchannelGrid(AMP::Mesh::Mesh::shared_ptr);   // Function to fill the subchannel data for all processors
    size_t N_subchannels;                                   // The total number of subchannels
    std::vector<double> d_x, d_y, d_z;                      // The x, y, z grid for the subchannel
    std::vector<bool> d_ownSubChannel;                      // Which subchannels do I own (multple procs my own a subchannel)
    std::vector<std::vector<int> > d_subchannelRanks;       // The processors that own each x-y point in the subchannel
    std::vector<std::vector<int> > d_subchannelRecv;        // The processors that are recieving data from each subchannel

    // Iterators over the mesh elemens of interest
    static AMP::Mesh::MeshIterator getSubchannelIterator(AMP::Mesh::Mesh::shared_ptr);
    AMP::Mesh::MeshIterator d_iterator1;
    AMP::Mesh::MeshIterator d_iterator2;

    // The list of local MeshElements in each subchannel
    std::vector<std::vector<AMP::Mesh::MeshElementID> >   d_elem;

    // Buffers to send/recv the data
    std::vector<MPI_Request> d_currRequests;
    std::vector<std::vector<double> > d_sendBuffer;

    int getSubchannelIndex( double x, double y );
};


} // Operator namespace
} // AMP namespace


