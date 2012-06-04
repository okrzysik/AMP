#ifndef included_AMP_MeshParameters
#define included_AMP_MeshParameters

#include "utils/AMP_MPI.h"
#include "utils/Database.h"


namespace AMP { 
namespace Mesh {



/**
 * \class Mesh
 * \brief A class used to pass the parameters for creating a mesh
 */
class MeshParameters
{
public:
    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a mesh manager object.
     */
    typedef boost::shared_ptr<MeshParameters>  shared_ptr;

    //! Empty constructor
    MeshParameters( );

    /**
     * \brief  Constructor to create the MeshParamaters from an AMP Database object.
     * \details  Constructor to create the MeshParamaters from an AMP Database object.
     * \param db    Input database for constructing a mesh
     */
    MeshParameters ( const boost::shared_ptr<AMP::Database>  db );

    /**
     * \brief       Set the comm for the mesh
     * \details     Set the communicator that will be used to construct the mesh
     * \param comm  The desired communicator
     */
    void setComm ( AMP::AMP_MPI comm );

    //!  Get the database for the mesh
    boost::shared_ptr<AMP::Database> getDatabase ( );

    //! Deconstructor
     ~MeshParameters ();

protected:

    //! A pointer to an AMP database containing the mesh info
    boost::shared_ptr<AMP::Database>  d_db;

    //! The desired communicator
    AMP::AMP_MPI comm;

    //! The maximum ghost cell width for the mesh
    int MAX_GCW_WIDTH;

//! See AMP::Mesh::Mesh for Mesh class
friend class Mesh;

};


} // Mesh namespace
} // AMP namespace

#endif

