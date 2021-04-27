#ifndef included_AMP_Map3to1to3
#define included_AMP_Map3to1to3

#include <map>
#include <vector>

#include "AMP/ampmesh/MeshIterator.h"
#include "AMP/operators/map/AsyncMapOperator.h"
#include "AMP/operators/map/AsyncMapOperatorParameters.h"


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
    explicit Map3to1to3( std::shared_ptr<const OperatorParameters> params );

    //!  Destructor
    virtual ~Map3to1to3();

    //! Return the name of the operator
    std::string type() const override { return "Map3to1to3"; }

    /** \brief   Set a frozen vector for results of the apply operation.
     * \details  Set a frozen vector for results of the apply operation.
     * \param result    The results vector
     */
    void setVector( AMP::LinearAlgebra::Vector::shared_ptr result ) override;

    /** \brief   Start a communicative apply operation.
     * \details  Start a communicative apply operation.
     */
    void applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                     AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /** \brief   Finish a communicative apply operation.
     * \details  Finish a communicative apply operation.
     */
    void applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                      AMP::LinearAlgebra::Vector::shared_ptr f ) override;


protected:
    /** \brief  Add an ordered pair to the set of interpolant values
     * \param map           The map to add the point to
     * \param [in] z        The ordinate
     * \param [in] val      The abscissa
     */
    void addTo1DMap( std::multimap<double, double> &map, double z, double val );

    /** \brief  Add an ordered pair to the set of interpolant values
     * \param map           The map to add the point to
     * \param [in] z        The ordinate
     * \param [in] val      The abscissa
     */
    void addTo1DMap( std::multimap<double, double> &map,
                     const std::vector<double> &z,
                     const std::vector<double> &val );

    /** \brief   A virtual method to construct the map from a vector
     * \details  This function constructs the map from a given vector.
     *    The inherited class must implement this function
     * \param [in] vec  The vector to be used to construct the map
     * \param [in] mesh The meshused to construct the map
     * \param [in] it   The iterator over the boundary used for the map
     */
    virtual std::multimap<double, double>
    buildMap( AMP::LinearAlgebra::Vector::const_shared_ptr vec,
              const AMP::Mesh::Mesh::shared_ptr mesh,
              const AMP::Mesh::MeshIterator &it );

    /** \brief  A virtual method to construct a vector from a map
     * \details  This function constructs a vector from the map.
     *    The inherited class must implement this function
     * \param [out] vec The vector to be used to construct the map
     * \param [in] mesh The mesh used to construct the map
     * \param [in] it   The iterator over the boundary used for the map
     * \param [in] map  The map containing all of the points
     */
    virtual void buildReturn( AMP::LinearAlgebra::Vector::shared_ptr vec,
                              const AMP::Mesh::Mesh::shared_ptr mesh,
                              const AMP::Mesh::MeshIterator &it,
                              const std::map<double, double> &map );

    //!  Iterators over the nodes on the boundary
    AMP::Mesh::MeshIterator d_srcIterator1;
    AMP::Mesh::MeshIterator d_srcIterator2;
    AMP::Mesh::MeshIterator d_dstIterator1;
    AMP::Mesh::MeshIterator d_dstIterator2;

private:
    //!  The place to put the mapped values
    AMP::LinearAlgebra::Vector::shared_ptr d_ResultVector;

    //!  A list of processors that own each mesh
    std::vector<bool> d_own_mesh1;
    std::vector<bool> d_own_mesh2;

    //!  The tag used for communication
    int d_commTag;

    //! structure used for communication
    struct comm_data {
        int N;
        double z;
        double sum;
        comm_data()
        {
            N   = 0;
            z   = 0.0;
            sum = 0.0;
        }
    };

    //!  The buffer used to perform the asynchronous communication
    std::vector<comm_data> d_SendBuf1;
    std::vector<comm_data> d_SendBuf2;

    // Function to unpack the recv buffer
    static void unpackBuffer( const std::vector<comm_data> &,
                              std::map<double, std::pair<int, double>> & );
};


} // namespace Operator
} // namespace AMP

#endif
