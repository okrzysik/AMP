#ifndef included_AMP_DOF_Manager
#define included_AMP_DOF_Manager

namespace AMP {
namespace Discretization {


/**
 * \class DOF_Manager
 * \brief A class used to provide DOF and vector creation routines
 *
 * \details  This class provides routines for calculating, accessing, and 
 *    using the degrees of freedom (DOF) per object.  It is also responsible 
 *    for creating vectors.
 */
class DOFManagerParameters
{
public:
    /**
     * \brief   Empty constructor for a DOF manager object
     * \details Empty constructor for a DOF manager object
     * \fn      DOFManagerParameters( )
     */
    DOFManagerParameters( );

    //! Deconstructor
     ~DOFManagerParameters ();


protected:

    //! Pointer to the underlying Mesh (may be NULL)
    boost::shared_ptr<AMP::Mesh::Mesh>  mesh;

    //! Pointer to the underlying VectorSpace (may be NULL)
    boost::shared_ptr<AMP::Discretization::VectorSpace>  vectorSpace;

};



} // Discretization namespace
} // AMP namespace

#endif

