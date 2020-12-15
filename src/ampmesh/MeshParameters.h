#ifndef included_AMP_MeshParameters
#define included_AMP_MeshParameters

#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"


namespace AMP {
namespace Mesh {


/**
 * \class Mesh
 * \brief A class used to pass the parameters for creating a mesh
 */
class MeshParameters
{
public:
    //! Empty constructor
    MeshParameters();

    /**
     * \brief  Constructor to create the MeshParamaters from an AMP Database object.
     * \details  Constructor to create the MeshParamaters from an AMP Database object.
     * \param db    Input database for constructing a mesh
     */
    explicit MeshParameters( const std::shared_ptr<AMP::Database> db );

    /**
     * \brief       Set the comm for the mesh
     * \details     Set the communicator that will be used to construct the mesh
     * \param comm  The desired communicator
     */
    void setComm( const AMP::AMP_MPI &comm );

    /**
     * \brief       Get the comm for the mesh
     * \details     Get the communicator for the mesh
     */
    const AMP::AMP_MPI &getComm() const { return comm; };

    //!  Get the database for the mesh
    std::shared_ptr<AMP::Database> getDatabase();

    //! Deconstructor
    virtual ~MeshParameters();

protected:
    //! A pointer to an AMP database containing the mesh info
    std::shared_ptr<AMP::Database> d_db;

    //! The desired communicator
    AMP::AMP_MPI comm;

    //! The maximum ghost cell width for the mesh
    int MAX_GCW_WIDTH;

    //! See AMP::Mesh::Mesh for Mesh class
    friend class Mesh;
};


} // namespace Mesh
} // namespace AMP

#endif
