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

    /**
     * \class  spMap
     * \brief  A class used to ...
     * \details  ...
     */
    class spMap
    {
    public:
        class type : public boost::shared_ptr<std::multimap<double,double> > 
        {
        public:
            bool isComputed;
            type ( ) : boost::shared_ptr<std::multimap<double,double> > ( ) {}
            type ( std::multimap<double,double> *a ) : boost::shared_ptr<std::multimap<double,double> > ( a ) {}
        };

        typedef std::multimap<double,double>::iterator iterator;

    protected:
        type   d_Map;

    public:
        void      insert ( std::pair<double,double> a ) { d_Map->insert ( a ); }
        iterator  find ( double a ) { return d_Map->find ( a ); }
        iterator  begin ()          { return d_Map->begin (); }
        iterator  end ()            { return d_Map->end (); }
        size_t    size ()           { return d_Map->size (); }
        iterator  upper_bound ( double a ) { return d_Map->upper_bound ( a ); }
        iterator  lower_bound ( double a ) { return d_Map->lower_bound ( a ); }
        void      clear () { d_Map->clear(); }
        void      erase ( iterator a , iterator b ) { d_Map->erase ( a , b ); }
        type     &data ()           { return d_Map; }
    };

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

    /** \brief   Continue the construction of the object requiring asynchronous calls. 
     * \details  Continue the construction of the object requiring asynchronous calls. 
     * \param[in] params  Input parameters
     */
    virtual bool continueAsynchronousConstruction ( const boost::shared_ptr < OperatorParameters > &params );

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

    /** \brief    
     */
    virtual AMP::LinearAlgebra::Variable::shared_ptr  getInputVariable( int varID = -1 ) { return d_MapVariable; }

    /** \brief   
     */
    virtual AMP::LinearAlgebra::Variable::shared_ptr  getOutputVariable() { return d_MapVariable; }

    /** \brief   
     */
    spMap::type  &getMapData () { return d_1DMultiMap.data(); }


protected:
    //!  An iterator over the nodes on the boundary
    AMP::Mesh::MeshIterator       d_BoundaryNodeIterator;

    //!  The variable of the appropriate vector to store the answer and to subset on the input
    AMP::LinearAlgebra::Variable::shared_ptr d_MapVariable;

    //!  The place to put the mapped values
    AMP::LinearAlgebra::Vector::shared_ptr d_ResultVector;

    //!  The storage for the interpolant
    spMap d_1DMultiMap;

    //!  A list of processors to exchange data with
    std::vector<char> d_SendToProc;

    //!  The buffer used to perform the asynchronous communication
    std::vector<double> d_SendBuf;

    //!  The communicator this map resides on
    AMP_MPI d_MapComm;

    //!  True if the map hasn't been "applied" yet
    bool d_FirstApply;

    //!  The tag used to send data
    int d_SendTag;

    //!  The tag used to recv data
    int d_RecvTag;

    /** \brief  This method will average nearby points in the interpolant
     * \param[in] tolerance  The distance across which values are assumed to be the same
     */
    virtual void smear ( double tolerance = 1.e-8 );

    /** \brief  Add an ordered pair to the set of interpolant values
     * \param [in] z   The ordinate
     * \param [in] val  The abscissa
     */
    void addTo1DMap ( double z , double val );

    /** \brief   A virtual method to construct the map from a vector
     * \details  This function constructs the map from a given vector.
     *    The inherited class must impliment this function
     * \param [in] p  The vector to be used to construct the map
     */
    virtual void buildMap ( const AMP::LinearAlgebra::Vector::shared_ptr p );

    /** \brief  A virtual method to construct a vector from a map
     * \details  This function constructs a vector from the map.
     *    The inherited class must impliment this function
     * \param [out] p  The vector to be constructed from the map
     */
    virtual void buildReturn ( AMP::LinearAlgebra::Vector::shared_ptr p );

};


} // Operator namespace
} // AMP namespace

#endif
