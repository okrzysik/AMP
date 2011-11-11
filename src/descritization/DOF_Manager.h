#ifndef included_AMP_DOF_Manager
#define included_AMP_DOF_Manager

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "descritization/DOF_ManagerParameters.h"
#include "vectors/Vector.h"
#include "matrices/Matrix.h"


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
class DOFManager
{
public:

    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a DOF manager object.
     */
    typedef boost::shared_ptr<AMP::Discretization::DOFManager>  shared_ptr;


    /** \brief Get the entry indices of nodal values given a mesh element
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     * \param[in]  obj      The element to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] ids      The entries in the vector associated with D.O.F.s on the nodes
     * \param[in]  which    Which D.O.F. to get.  If not specified, return all D.O.F.s
     */
    virtual void getDOFs ( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &ids , unsigned int which = static_cast<unsigned int>(-1) ) const;


    /** \brief  The first D.O.F. on this core
     * \return The first D.O.F. on this core
     */
    virtual size_t  beginDOF ( );


    /** \brief  One past the last D.O.F. on this core
     * \return One past the last D.O.F. on this core
     */
    virtual size_t  endDOF ( );


    /** \brief  The local number of D.O.F 
     * \return  The local number of D.O.F 
     */
    virtual size_t  numLocalDOF ( );


    /** \brief  The global number of D.O.F 
     * \return  The global number of D.O.F 
     */
    virtual size_t  numGlobalDOF ( );


    /**
     * \brief Create a new AMP vector
     * \details  This function creates a new AMP vector for the given variable, using the current DOF properties.
     * \param variable  Variable that will be used to create the vector
     */
    virtual AMP::LinearAlgebra::Vector::shared_ptr   createVector ( AMP::LinearAlgebra::Variable::shared_ptr variable );


    /**
     * \brief Create a new AMP matrix
     * \details  This function creates a new AMP matrix for the given variable, using the current DOF properties.
     * \param operand  Variable that will be used to create the matrix
     * \param result   Variable that will be used to create the matrix
     */
    virtual  AMP::LinearAlgebra::Matrix::shared_ptr   createMatrix ( AMP::LinearAlgebra::Variable::shared_ptr operand , AMP::LinearAlgebra::Variable::shared_ptr result );

 

protected:

    //!  Empty constructor for a DOF manager object
    DOFManager ( ) {};

    //! The DOF manager parameters
    const DOFManagerParameters::shared_ptr params;

};



} // Discretization namespace
} // AMP namespace

#endif

