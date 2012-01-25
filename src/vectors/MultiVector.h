#ifndef included_AMP_MultiVector
#define included_AMP_MultiVector

#include "Vector.h"
#include "MultiVariable.h"
#include "VectorEngine.h"
#include "DataChangePassThrough.h"

namespace AMP {
namespace LinearAlgebra {

  /// \brief Typedef for the parameters to create a MultiVector
  typedef VectorParameters MultiVectorParameters;

  /** \brief A collection of AMP Vectors that appear as one vector
    * \details
       Given a set of vectors, they can be collected into a single vector.  This class
       accomplishes this task.
       */
  class MultiVector : public Vector , public VectorEngine , public DataChangePassThrough
  {
    protected:
      /// \brief  The communicator this multivector is on
      AMP_MPI                            d_Comm;
      /// \brief  Indicates if the multivector created the communicator
      bool                               d_CommCreated;

      /// \brief The list of AMP Vectors that comprise this Vector
      std::vector<Vector::shared_ptr>    d_vVectors;
      /// \brief The global lengths of each vector
      std::vector<size_t>                d_vGlobalOffsets;
      /// \brief The local length of each vector
      std::vector<size_t>                d_vLocalOffsets;

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
            Vector::shared_ptr &getVector (       VectorOperations &vec , size_t which ) const;

      /** \brief A method that will translate an array of global ids relative to the multivector
        * into an array of arrays of global ids relative to the component vectors
        * \param[in] num            The length of the id array indices and double array vals
        * \param[in] indices        The indices of the values relative to the multivector
        * \param[in] vals           Values associated somehow with the indices
        * \param[out] out_indices   An array of arrays of mapped indices relative to constituent vectors
        * \param[out] out_vals      The values partitioned according to out_indices
        * \param[in] offset_vec     The lengths of the constituent vectors
        * \param[out] remap         If not null, this is a list of where in the indices array 
        *                           an entry comes from
        * \details  This is the method that handles all of the bookkeeping associated with
        * treating a set of vectors of arbitrary length as a single vector
        */
      void  partitionValues ( int num , 
                              int *indices , 
                              const double *vals ,
                              std::vector<std::vector<int> >  &out_indices , 
                              std::vector<std::vector<double > >  &out_vals ,
                              const std::vector<size_t>  &offset_vec , 
                              std::vector< std::vector<size_t> >  *remap = 0 ) const ;

      /** Constructor:  create a MultiVector with a particular variable
        * \param[in]  names  The vector to create the MultiVector from
        */
      MultiVector ( Variable::shared_ptr names );

      virtual void dataChanged ();

    public:
      /** Iterator typedef
        */
      typedef std::vector<Vector::shared_ptr>::iterator  vector_iterator;

      /** Iterator typedef
        */
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
        */
      static Vector::shared_ptr  create ( Variable::shared_ptr name , AMP_MPI  comm );

      /** \brief Create a new multivector in parallel
        * \param[in] name  Name of the new vector
        * \param[in] comm  Communicator to build the MultiVector on
        */
      static Vector::shared_ptr  create ( const std::string &name , AMP_MPI comm );

      /** \brief Create a multivector view of a vector
        * \param[in] vec  The vector to view
        * \param[in] comm  Communicator to create the MultiVector on
        * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
        * and vec is added to it.  If vec is not a parallel vector (such as a SimpleVector), comm
        * must be specified.
        */
      static Vector::shared_ptr  view ( Vector::shared_ptr &vec , AMP_MPI comm = AMP_MPI(AMP_COMM_NULL) );

      /** \brief Encapsulate a vector in a MultiVector
        * \param[in] vec  The vector to view
        * \param[in] comm  Communicator to create the MultiVector on
        * \details  If vec is a MultiVector, it is returned.  Otherwise, a MultiVector is created
        * and vec is added to it.  If vec is not a parallel vector (such as a SimpleVector), comm
        * must be specified.
        */
      static Vector::shared_ptr  encapsulate ( Vector::shared_ptr &vec , AMP_MPI comm = AMP_MPI(AMP_COMM_NULL) );

      virtual void addVector ( Vector::shared_ptr );

      /** \brief  Remove a vector from a MultiVector
        * \param[in]  vec  The Vector to remove from the MultiVector
        * \return  The next Vector in the MultiVector
        */
      virtual vector_iterator eraseVector ( vector_iterator vec );

      /** \brief  Return the i-th Vector in this MultiVector
        * \param[in]  i  Which vector to return
        * \return  The i-th Vector
        */
      virtual Vector::shared_ptr  getVector ( size_t i );

      /** \brief  Obtain the number of Vector objects in this MultiVector
        * \return  The number of Vector objects in this MultiVector
        */
      size_t  getNumberOfSubvectors ();

      /** \brief  Return the beginning iterator of the i-th Vector
        * \param[in] i  Which Vector to return the beginning iterator from
        * \return The iterator
        */
      virtual Vector::iterator    getIterator ( size_t i );

      /** \brief  Return the beginning iterator of the i-th Vector
        * \param[in] i  Which Vector to return the beginning iterator from
        * \return The iterator
        */
      virtual Vector::const_iterator    getIterator ( size_t i ) const;

      /** \brief  Destructor
        */
      virtual ~MultiVector();

      virtual Vector::iterator    begin();
      virtual Vector::iterator    end();
      virtual Vector::const_iterator    begin() const;
      virtual Vector::const_iterator    end() const;

      virtual void        selectInto ( const VectorSelector & , Vector::shared_ptr );

      virtual void        dumpOwnedData ( std::ostream &out , size_t GIDoffset=0 , size_t LIDoffset = 0 ) const;
      virtual void        dumpGhostedData ( std::ostream &out , size_t offset=0 ) const;
      virtual std::string type() const;

      virtual AMP_MPI  getComm () const;

      virtual size_t numberOfDataBlocks () const;

      virtual size_t sizeOfDataBlock ( size_t i ) const;

      virtual Vector::shared_ptr  subsetVectorForVariable ( const Variable::shared_ptr &name );
      virtual Vector::shared_ptr cloneVector(const Variable::shared_ptr name) const;
//      virtual void copyVector(const Vector &src_vec);
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
      virtual void setValuesByLocalID ( int num , int *indices , const double *vals );
      virtual void setLocalValuesByGlobalID ( int num , int *indices , const double *vals );
      virtual void setValuesByGlobalID ( int num , int *indices , const double *vals );
      virtual void addValuesByLocalID ( int num , int *indices , const double *vals );
      virtual void addLocalValuesByGlobalID ( int num , int *indices , const double *vals );
      virtual void addValuesByGlobalID ( int num , int *indices , const double *vals );
      virtual void getValuesByGlobalID ( int numVals , int *ndx , double *vals ) const;
      virtual void getLocalValuesByGlobalID ( int numVals , int *ndx , double *vals ) const;
      virtual void getValuesByLocalID ( int numVals , int *ndx , double *vals ) const;
      virtual void makeConsistent ( ScatterType  t );
      virtual void assemble();
      virtual double L1Norm(void) const;
      virtual double L2Norm(void) const;
      virtual double maxNorm(void) const;
      virtual double dot(const VectorOperations &x) const;
      virtual size_t getLocalSize() const;
      virtual size_t getGlobalSize() const;
      virtual size_t getGhostSize() const;
      virtual void   putRawData ( double * );


      ///Vector engine functions
      virtual VectorEngine::BufferPtr  getNewBuffer();
      virtual bool                     sameEngine ( VectorEngine &rhs ) const;
      virtual VectorEngine::shared_ptr cloneEngine ( VectorEngine::BufferPtr  ) const;
      virtual void                     swapEngines ( VectorEngine::shared_ptr p );
      virtual const void              *getDataBlock ( size_t i ) const;
      virtual void                    *getDataBlock ( size_t i );

      virtual void   copyOutRawData ( double ** out );
  };

}
}

#include "MultiVector.inline.h"


#endif


