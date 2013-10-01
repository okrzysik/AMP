#ifndef included_AMP_MultiVector
#define included_AMP_MultiVector

#include "Vector.h"
#include "MultiVariable.h"
#include "VectorEngine.h"
#include "DataChangePassThrough.h"

namespace AMP {
namespace LinearAlgebra {


/** \brief A collection of AMP Vectors that appear as one vector
  * \details
  *    Given a set of vectors, they can be collected into a singlevector.  This class
  *    accomplishes this task.
  */
class MultiVector : public Vector , public VectorEngine , public DataChangePassThrough
{
public:
    //!  Iterator typedef
    typedef std::vector<Vector::shared_ptr>::iterator  vector_iterator;

    //!  Iterator typedef
    typedef std::vector<Vector::shared_ptr>::const_iterator  vector_const_iterator;

    /** \brief Return the first vector in the MultiVector
      * \return The first vector
      */
    vector_iterator beginVector();

    /** \brief Return one past the last vector in the MultiVector
      * \return One past the last vector
      */
    vector_iterator endVector();

    /** \brief Determine if a Vector is a constituent
      * \param[in]  p  The vector to look for
      * \return True if the pointer p is in the multivector
      */
    bool  containsPointer ( const Vector::shared_ptr p ) const;

    /** \brief Create a new multivector in parallel
      * \param[in] name  Variable describing the new vector
      * \param[in] comm  Communicator to build the MultiVector on
      * \param[in] vec   Optional list of vectors in the MultiVector
      */
    static boost::shared_ptr<MultiVector>  create ( Variable::shared_ptr name, AMP_MPI comm, 
        const std::vector<Vector::shared_ptr>& vecs=std::vector<Vector::shared_ptr>() );

    /** \brief Create a new multivector in parallel
      * \param[in] name  Name of the new vector
      * \param[in] comm  Communicator to build the MultiVector on
      * \param[in] vec   Optional list of vectors in the MultiVector
      */
    static boost::shared_ptr<MultiVector>  create ( const std::string &name, AMP_MPI comm, 
        const std::vector<Vector::shared_ptr>& vecs=std::vector<Vector::shared_ptr>() );

    /** \brief Create a new multivector in parallel
      * \param[in] name  Variable describing the new vector
      * \param[in] comm  Communicator to build the MultiVector on
      * \param[in] vec   Optional list of vectors in the MultiVector
      */
    static boost::shared_ptr<const MultiVector>  const_create ( Variable::shared_ptr name, AMP_MPI comm, 
        const std::vector<Vector::const_shared_ptr>& vecs=std::vector<Vector::const_shared_ptr>() );

    /** \brief Create a new multivector in parallel
      * \param[in] name  Name of the new vector
      * \param[in] comm  Communicator to build the MultiVector on
      * \param[in] vec   Optional list of vectors in the MultiVector
      */
    static boost::shared_ptr<const MultiVector>  const_create ( const std::string &name, AMP_MPI comm, 
        const std::vector<Vector::const_shared_ptr>& vecs=std::vector<Vector::const_shared_ptr>() );

    /** \brief Create a multivector view of a vector
      * \param[in] vec  The vector to view
      * \param[in] comm  Communicator to create the MultiVector on
      * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
      * and vec is added to it.  If vec is not a parallel vector (such as a SimpleVector), comm
      * must be specified.
      */
    static boost::shared_ptr<MultiVector>   view ( Vector::shared_ptr &vec , AMP_MPI comm = AMP_MPI(AMP_COMM_NULL) );

    /** \brief Encapsulate a vector in a MultiVector
      * \param[in] vec  The vector to view
      * \param[in] comm  Communicator to create the MultiVector on
      * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
      * and vec is added to it.  If vec is not a parallel vector (such as a SimpleVector), comm
      * must be specified.
      */
    static boost::shared_ptr<MultiVector>   encapsulate ( Vector::shared_ptr &vec , AMP_MPI comm = AMP_MPI(AMP_COMM_NULL) );

    /** \brief Replace a vector in a MultiVector
      * \details  This function will replace a given vector in the multivector with a different
      * vector.  The original and new vectors must share the same DOFManager. 
      * This is a purly local operation and does not require communication, but should be called 
      * by all processors that own the vectors.
      * \param[in] oldVec  The original vector to replace
      * \param[in] newVec  The new vector to use
      */
    virtual void replaceSubVector(Vector::shared_ptr oldVec, Vector::shared_ptr newVec) ;

    /** \brief  Add a vector to a MultiVector.
      *   Note: this is a collective operation, vec may be NULL on some processors
      * \param[in]  vec  The Vector to add to the MultiVector
      */
    virtual void addVector ( Vector::shared_ptr vec );

    /** \brief  Add vector(s) to a MultiVector.
      *   Note: This is a collective operation, vec may be different sizes on different processors
      * \param[in]  vec  The Vector to add to the MultiVector
      */
    virtual void addVector ( std::vector<Vector::shared_ptr> vec );

    /** \brief  Remove a vector from a MultiVector
      * \param[in]  vec  The Vector to remove from the MultiVector
      */
    virtual void eraseVector ( Vector::shared_ptr vec );

    /** \brief  Return the i-th Vector in this MultiVector
      * \param[in]  i  Which vector to return
      * \return  The i-th Vector
      */
    virtual Vector::shared_ptr  getVector ( size_t i );

    /** \brief  Return the i-th Vector in this MultiVector
      * \param[in]  i  Which vector to return
      * \return  The i-th Vector
      */
    virtual Vector::const_shared_ptr  getVector ( size_t i ) const;

    /** \brief  Obtain the number of Vector objects in this MultiVector
      * \return  The number of Vector objects in this MultiVector
      */
    size_t  getNumberOfSubvectors () const;

    //!  Destructor
    virtual ~MultiVector();

    // Vector functions
    using Vector::cloneVector;

    virtual void      dumpOwnedData ( std::ostream &out , size_t GIDoffset=0 , size_t LIDoffset = 0 ) const;
    virtual void      dumpGhostedData ( std::ostream &out , size_t offset=0 ) const;
    virtual std::string type() const;

    virtual AMP_MPI  getComm () const;

    virtual size_t numberOfDataBlocks () const;

    virtual size_t sizeOfDataBlock ( size_t i ) const;

    virtual Vector::shared_ptr  subsetVectorForVariable ( const Variable::shared_ptr &name );
    virtual Vector::const_shared_ptr  constSubsetVectorForVariable ( const Variable::shared_ptr &name ) const;
    virtual Vector::shared_ptr cloneVector(const Variable::shared_ptr name) const;
    virtual void copyVector(const Vector::const_shared_ptr &src_vec);
    virtual void swapVectors(Vector &other);
    virtual void aliasVector(Vector &other);
    virtual void setToScalar(double alpha);
    virtual void scale(double alpha, const VectorOperations &x);
    virtual void scale(double alpha);
    virtual void add(const VectorOperations &x, const VectorOperations &y);
    virtual void subtract(const VectorOperations &x, const VectorOperations &y);
    virtual void multiply( const VectorOperations &x, const VectorOperations &y);
    virtual void divide( const VectorOperations &x, const VectorOperations &y);
    virtual void reciprocal(const VectorOperations &x);
    virtual void linearSum(double alpha, const VectorOperations &x,
          double beta, const VectorOperations &y);
    virtual void axpy(double alpha, const VectorOperations &x, const VectorOperations &y);
    virtual void axpby(double alpha, double beta, const VectorOperations &x);
    virtual void abs(const VectorOperations &x);
    virtual double min(void) const;
    virtual double max(void) const;
    virtual void setRandomValues(void);
    virtual void setValuesByLocalID ( int num , size_t *indices , const double *vals );
    virtual void setLocalValuesByGlobalID ( int num , size_t *indices , const double *vals );
    virtual void setGhostValuesByGlobalID ( int num , size_t *indices , const double *vals );
    virtual void setValuesByGlobalID ( int num , size_t *indices , const double *vals );
    virtual void addValuesByLocalID ( int num , size_t *indices , const double *vals );
    virtual void addLocalValuesByGlobalID ( int num , size_t *indices , const double *vals );
    virtual void addValuesByGlobalID ( int num , size_t *indices , const double *vals );
    virtual void getValuesByGlobalID ( int numVals , size_t *ndx , double *vals ) const;
    virtual void getLocalValuesByGlobalID ( int numVals , size_t *ndx , double *vals ) const;
    virtual void getGhostValuesByGlobalID ( int numVals , size_t *ndx , double *vals ) const;
    virtual void getValuesByLocalID ( int numVals , size_t *ndx , double *vals ) const;
    virtual void makeConsistent ( ScatterType  t );
    virtual UpdateState  getUpdateStatus() const;
    virtual void  setUpdateStatus( UpdateState state );
    virtual void assemble();
    virtual double L1Norm(void) const;
    virtual double L2Norm(void) const;
    virtual double maxNorm(void) const;
    using Vector::dot;
    virtual double dot(const VectorOperations &x) const;
    virtual size_t getLocalSize() const;
    virtual size_t getGlobalSize() const;
    virtual size_t getGhostSize() const;
    virtual void   putRawData ( double * );

    // Vector engine functions
    virtual boost::shared_ptr<std::vector<double> >  getNewBuffer();
    virtual bool               sameEngine ( VectorEngine &rhs ) const;
    virtual VectorEngine::shared_ptr cloneEngine ( boost::shared_ptr<std::vector<double> > ) const;
    virtual void               swapEngines ( VectorEngine::shared_ptr p );
    virtual const void          *getDataBlock ( size_t i ) const;
    virtual void              *getDataBlock ( size_t i );

    virtual void   copyOutRawData ( double ** out );


protected:

    virtual void selectInto ( const VectorSelector & , Vector::shared_ptr );
    virtual void constSelectInto ( const VectorSelector &criterion , Vector::shared_ptr vector ) const;

    //! The communicator this multivector is on
    AMP_MPI                    d_Comm;
    //! Indicates if the multivector created the communicator
    bool                     d_CommCreated;

    //! The list of AMP Vectors that comprise this Vector
    std::vector<Vector::shared_ptr>    d_vVectors;

    virtual void *getRawDataBlockAsVoid ( size_t );
    virtual const void *getRawDataBlockAsVoid ( size_t ) const;

    /** \brief A convenience method for extracting vectors from a base class
      * \param[in]  vec  The vector to extract a vector from
      * \param[in]  which  Which vector to get
      * \return     The extracted vector
      */
    const Vector::shared_ptr &getVector ( const VectorOperations &vec , size_t which ) const;

    /** \brief A convenience method for extracting vectors from a base class
      * \param[in]  vec  The vector to extract a vector from
      * \param[in]  which  Which vector to get
      * \return     The extracted vector
      */
    Vector::shared_ptr &getVector (     VectorOperations &vec , size_t which ) const;

    /** Constructor:  create a MultiVector with a particular variable
      * \param[in]  names  The vector to create the MultiVector from
      */
    MultiVector ( const std::string& name );

    virtual void dataChanged ();

    /** A method that will translate an array of global ids relative to the multivector
      * into an array of arrays of global ids relative to the component vectors
      * \param[in] num            The number of DOFs that need to be mapped
      * \param[in] indices        The indices of the values relative to the multivector
      * \param[in] vals           Values associated somehow with the indices
      * \param[out] out_indices   An array of arrays of mapped indices relative to constituent vectors
      * \param[out] out_vals      The values partitioned according to out_indices
      * \param[out] remap         If not null, this is a list of where in the indices array an entry comes from
      */
    void  partitionGlobalValues ( const int num, 
                            const size_t *indices, 
                            const double *vals,
                            std::vector<std::vector<size_t> >  &out_indices, 
                            std::vector<std::vector<double> >  &out_vals,
                            std::vector<std::vector<int> > *remap=NULL ) const;

    /** A method that will translate an array of local ids relative to the multivector
      * into an array of arrays of local ids relative to the component vectors
      * \param[in] num            The number of DOFs that need to be mapped
      * \param[in] indices        The indices of the values relative to the multivector
      * \param[in] vals           Values associated somehow with the indices
      * \param[out] out_indices   An array of arrays of mapped indices relative to constituent vectors
      * \param[out] out_vals      The values partitioned according to out_indices
      * \param[out] remap         If not null, this is a list of where in the indices array an entry comes from
      */
    void  partitionLocalValues ( const int num, 
                            const size_t *indices, 
                            const double *vals,
                            std::vector<std::vector<size_t> >  &out_indices, 
                            std::vector<std::vector<double> >  &out_vals,
                            std::vector<std::vector<int> > *remap=NULL ) const;

};


}
}

#include "MultiVector.inline.h"


#endif


