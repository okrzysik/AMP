#ifndef included_AMP_DualVector
#define included_AMP_DualVector

#include "Vector.h"
#include "VectorEngine.h"
#include "DualVariable.h"
#include "DataChangePassThrough.h"

namespace AMP {
namespace LinearAlgebra {

  /**
    * \class DualVectorParameters
    * \brief Parameters used to create a DualVector
    */
  class DualVectorParameters : public VectorParameters {};

  /**
    * \brief A vector composed of two vectors
    * \details  Let \f$ v = [ v_1 \  v_2 ]^T \f$.  In this instance,
    * \f$v\f$ would be a DualVector made up of two vectors \f$v_1\f$ and
    * \f$v_2\f$.  The constituents of a DualVector do not need to be
    * contiguous in memory.  In fact, the constituents need not be
    * contiguous themselves, they must only present the contiguous interface
    * of Vector.
    *
    * In order to present a contiguous interface, the indices in the vector
    * are given as \f$ v = [ v_1 \  v_2 ]^T \f$.  In other words, the 
    * entries in \f$v\f$ which point to \f$v_1\f$ are numbered as \f$v_1\f$.
    * The entries associated with \f$v_2\f$ are offset by the number of
    * values in \f$v_1\f$.
    *
    * For instance, let \f$ v_1 = [ 1\ 2\ 3\ 4\ 5 ]^T\f$ and 
    * \f$ v_2 = [ 10\ 11\ 12 ]^T \f$.  Then, 
    * \f$ v = [ 1\ 2\ 3\ 4\ 5\ 10\ 11\ 12 ]^T \f$.
    * 
    * DualVector is also a VectorEngine allowing for views of DualVector through
    * a ManagedVector wrapper.
    *
    * This is a class that is used for composing vectors from other vectors.
    * Only use this class to create the composition.  Once composed, there is
    * no reason to cast the vector to this type.  For instance, the only
    * valid use of DualVector in code would be
    * \code
      Vector::shared_ptr v1;
      Vector::shared_ptr v2;
         .
         .
         .
      Vector::shared_ptr v = DualVector::create ( v1 , v2, "New Dual Vector" );
      \endcode
    */
  class DualVector : public Vector , 
                     public VectorEngine , 
                     public DataChangePassThrough
  {
    protected:
      /// \cond UNDOCUMENTED
      /**
        * \brief The first vector in the DualVector
        */
      Vector::shared_ptr     d_pVector1;

      /**
        * \brief The second vector in the DualVector
        */
      Vector::shared_ptr     d_pVector2;

      void *getRawDataBlockAsVoid ( size_t blockNum );
      const void *getRawDataBlockAsVoid ( size_t blockNum) const;

      /**
        * \brief This will return the first vector in rhs
        * \param rhs  A VectorOperations that is a DualVector
        * \details  This convenience method will attempt to cast rhs to a 
        * DualVector.  If this cast succeeds, then it will return the first
        * vector in rhs.
        */
      static const Vector::shared_ptr &getVector1 ( const VectorOperations &rhs );
      
      /**
        * \brief This will return the second vector in rhs
        * \param rhs  A VectorOperations that is a DualVector
        * \details  This convenience method will attempt to cast rhs to a 
        * DualVector.  If this cast succeeds, then it will return the second
        * vector in rhs.
        */
      static const Vector::shared_ptr &getVector2 ( const VectorOperations &rhs );

      /**
        * \brief This will return the first vector in rhs
        * \param rhs  A VectorOperations that is a DualVector
        * \details  This convenience method will attempt to cast rhs to a 
        * DualVector.  If this cast succeeds, then it will return the first
        * vector in rhs.
        */
      static Vector::shared_ptr &getVector1 ( VectorOperations &rhs );
      
      /**
        * \brief This will return the second vector in rhs
        * \param rhs  A VectorOperations that is a DualVector
        * \details  This convenience method will attempt to cast rhs to a 
        * DualVector.  If this cast succeeds, then it will return the second
        * vector in rhs.
        */
      static Vector::shared_ptr &getVector2 ( VectorOperations &rhs );

      /**
        * \brief Given a set of indices with global ids w.r.t. this DualVector,
        * this function will return indices for each constituent vector w.r.t.
        * each vector.  It also partitions the values among the constitutents
        * \param num  Number of incoming indices and values
        * \param indices  Indices given w.r.t. this DualVector
        * \param vals  The values to be partitioned
        * \param start2 The number of elements in the first vector, either
        * in a local sense or a global sense
        * \param Ndx1  The indices destined for the first vector
        * \param Ndx2  The indices destined for the second vector
        * \param Vals1 The values destined for the first vector
        * \param Vals2 The values destined for the second vector
        */
      void  partitionValues ( int num , int *indices , const double *vals , int start2 ,
                    std::vector<int> &Ndx1 , std::vector<int> &Ndx2 , 
                    std::vector<double> &Vals1 , std::vector<double> &Vals2 ) const;

      DualVector ( Vector::shared_ptr v1 , Vector::shared_ptr v2 , Variable::shared_ptr name );
      /// \endcond

    public:

      /**
        * \brief DualVector factory
        * \param[in] v1  First vector
        * \param[in] v2  Second vector
        * \param[in] name Name of the new Variable to associate with this vector
        * \return A shared pointer to a Vector
        */
      static Vector::shared_ptr  create ( Vector::shared_ptr v1 , Vector::shared_ptr v2 , const std::string &name );

      /**
        * \brief DualVector factory
        * \param[in] v1  First vector
        * \param[in] v2  Second vector
        * \param[in] name Variable to associate with this vector
        * \return A shared pointer to a Vector
        */
      static Vector::shared_ptr  create ( Vector::shared_ptr v1 , Vector::shared_ptr v2 , Variable::shared_ptr name );


      /// \cond UNDOCUMENTED
      virtual std::string type() const { return "Dual Vector"; }
      virtual void  selectInto ( const VectorSelector & , Vector::shared_ptr );
      virtual ~DualVector(); 
      virtual AMP_MPI  getComm () const;
      virtual size_t numberOfDataBlocks () const;
      virtual size_t sizeOfDataBlock ( size_t i ) const;
      virtual void * getDataBlock ( size_t i );
      virtual const void * getDataBlock ( size_t i ) const;
      virtual Vector::shared_ptr  subsetVectorForVariable ( const Variable::shared_ptr &name );

      virtual Vector::shared_ptr cloneVector(const Variable::shared_ptr name) const;
      virtual void copyVector(const Vector &src_vec);
      virtual void swapVectors(Vector &other);
      virtual void aliasVector(Vector &other);
      virtual size_t getLocalSize() const;
      virtual size_t getGlobalSize() const;
      virtual size_t getGhostSize() const;
      virtual void   putRawData ( double * );
      virtual void   copyOutRawData ( double **out );
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
      virtual double L1Norm(void) const;
      virtual double L2Norm(void) const;
      virtual double maxNorm(void) const;
      virtual double dot(const VectorOperations &x) const;

      virtual void setValuesByLocalID ( int , int * , const double * );
      virtual void setLocalValuesByGlobalID ( int , int * , const double * );
      virtual void setValuesByGlobalID ( int , int * , const double * );
      virtual void addValuesByLocalID ( int , int * , const double * );
      virtual void addLocalValuesByGlobalID ( int , int * , const double * );
      virtual void addValuesByGlobalID ( int , int * , const double * );
      virtual void getValuesByGlobalID ( int , int * , double * ) const;
      virtual void getLocalValuesByGlobalID ( int , int * , double * ) const;
      virtual void getValuesByLocalID ( int , int * , double * ) const;

      virtual void makeConsistent ( ScatterType  t );
      virtual void assemble();


      virtual BufferPtr  getNewBuffer ();
      virtual bool       sameEngine ( VectorEngine &p ) const;
      virtual void       swapEngines ( VectorEngine::shared_ptr p );
      virtual VectorEngine::shared_ptr  cloneEngine ( BufferPtr ) const;
      /// \endcond
  };

}
}

#include "DualVector.inline.h"

#endif


