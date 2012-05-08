#ifndef  included_AMP_Map3to1to3
#define  included_AMP_Map3to1to3

#include <map>
#include <vector>

#include "operators/map/AsyncMapOperator.h"
#include "operators/map/AsyncMapOperatorParameters.h"
#include "ampmesh/MeshIterator.h"


namespace AMP {
namespace Operator {

  /** \brief  For problems with a great deal of symmetry, it is possible
    * to move data from one mesh to another by generating a 1D appoximation
    * using one mesh and interpolating the results on another
    */


/**
 * \class  Map3to1to3
 * \brief  A class used to reduce a 3D problem to 1D, transfer the solution, and map back to 3D
 * \details  For problems with a great deal of symmetry, it is possible
 *   to move data from one mesh to another by generating a 1D appoximation
 *   using one mesh and interpolating the results on another.  This class 
 *   is a base class to manage the data transfer from 3D to 1D back to 3D.  
 *   An inherited map must impliment buildMap and buildReturn to manage how
 *   the data is mapped between 3D and 1D.
 */
class Map3to1to3 : public AsyncMapOperator
{
public:

    /** \brief   Standard constructor
     * \param[in] params  Input parameters
     */
    Map3to1to3 ( const boost::shared_ptr<OperatorParameters> & params );

    //!  Destructor
    virtual ~Map3to1to3 ();

    /** \brief   Set a frozen vector for results of the apply operation. 
     * \details  Set a frozen vector for results of the apply operation. 
     * \param result    The results vector
     */
    virtual void  setVector ( AMP::LinearAlgebra::Vector::shared_ptr &result );

    /** \brief   Start a communicative apply operation. 
     * \details  Start a communicative apply operation. 
     */
    virtual void applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const  AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r,
        const double a = -1.0, const double b = 1.0);

    /** \brief   Finish a communicative apply operation. 
     * \details  Finish a communicative apply operation. 
     */
    virtual void applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const  AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r,
        const double a = -1.0, const double b = 1.0);


protected:
    /** \brief  This method will average nearby points in the interpolant
     * \param map           The map to smear
     * \param [in] tolerance  The distance across which values are assumed to be the same
     */
    virtual void smear( std::multimap<double,double> &map, double tolerance=1.e-8 );

    /** \brief  Add an ordered pair to the set of interpolant values
     * \param map           The map to add the point to
     * \param [in] z        The ordinate
     * \param [in] val      The abscissa
     */
    void addTo1DMap( std::multimap<double,double> &map, double z , double val );

    /** \brief   A virtual method to construct the map from a vector
     * \details  This function constructs the map from a given vector.
     *    The inherited class must impliment this function
     * \param [in] vec  The vector to be used to construct the map
     * \param [in] mesh The meshused to construct the map
     * \param [in] it   The iterator over the boundary used for the map
     */
    virtual std::multimap<double,double>  buildMap( const AMP::LinearAlgebra::Vector::shared_ptr vec, 
        const AMP::Mesh::Mesh::shared_ptr mesh, const AMP::Mesh::MeshIterator &it );

    /** \brief  A virtual method to construct a vector from a map
     * \details  This function constructs a vector from the map.
     *    The inherited class must impliment this function
     * \param [out] vec The vector to be used to construct the map
     * \param [in] mesh The meshused to construct the map
     * \param [in] it   The iterator over the boundary used for the map
     * \param [in] map  The map containing all of the points
     */
    virtual void buildReturn( AMP::LinearAlgebra::Vector::shared_ptr vec, const AMP::Mesh::Mesh::shared_ptr mesh, 
        const AMP::Mesh::MeshIterator &it, const std::multimap<double,double> &map );

private:
    //!  Iterators over the nodes on the boundary
    AMP::Mesh::MeshIterator  d_iterator1;
    AMP::Mesh::MeshIterator  d_iterator2;

    //!  The place to put the mapped values
    AMP::LinearAlgebra::Vector::shared_ptr d_ResultVector;

    //!  A list of processors that own each mesh
    std::vector<bool> d_own_mesh1;
    std::vector<bool> d_own_mesh2;

    //!  The buffer used to perform the asynchronous communication
    std::vector<double> d_SendBuf1;
    std::vector<double> d_SendBuf2;

    //!  True if the map hasn't been "applied" yet
    bool d_FirstApply;

    //!  The tag used for communication
    int d_commTag;

};


} // Operator namespace
} // AMP namespace

#endif
