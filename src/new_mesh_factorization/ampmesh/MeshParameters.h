#ifndef included_AMP_MeshParameters
#define included_AMP_MeshParameters


namespace AMP { 
namespace Mesh {


//! Enumeration for basic mesh-based quantities
enum ElementType { Vertex=0, Edge=1, Face=2, Volume=3 };


/**
 * \class Mesh
 * \brief A class used to pass the parameters for creating a mesh
 */
class MeshParamters
{
public:
    /**
     * \brief  Constructor to create the MeshParamaters from an AMP Database object.
     * \details  Constructor to create the MeshParamaters from an AMP Database object.
     * \fn          Mesh ( const boost::shared_ptr<AMP::Database>  & db )
     * \param db    Input database for constructing a mesh
     */
    MeshParameters ( const boost::shared_ptr<AMP::Database>  & db );

    //! Deconstructor
     ~MeshParameters ();
protected:

    //! A pointer to an AMP database containing the mesh info
    boost::shared_ptr<AMP::Database>  d_db;

    //! The maximum ghost cell width for the mesh
    const int MAX_GCW_WIDTH;
};


} // Mesh namespace
} // AMP namespace

#endif

