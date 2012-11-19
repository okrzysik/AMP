#ifndef included_AMP_ManagedVector
#define included_AMP_ManagedVector


#include <vector>
#include <stdexcept>
#include "Vector.h"
#include "DataChangeFirer.h"
#include "MultiVector.h"


namespace AMP {
namespace LinearAlgebra {


/**
  \brief Data necessary to create a managed vector
*/
class ManagedVectorParameters : public VectorParameters
{
protected:
    //!  Copy constructor is protected to prevent unintended copies
    ManagedVectorParameters ( const ManagedVectorParameters & );

public:
    //! Constructor
    ManagedVectorParameters ();

    //! The VectorEngine to use with the managed vector
    VectorEngine::shared_ptr    d_Engine;

    //! Indicates whether the engine should be used as is or cloned
    bool            d_CloneEngine;

    //! Buffer to use for the managed vector
    VectorEngine::BufferPtr     d_Buffer;
};


/**
   \brief Class used to control data and kernels of various vector libraries
   \details  A ManagedVector will take an engine and create a buffer, if 
   necessary.  

   A ManagedVector has two pointers: data and engine.  If the data pointer
   is null, then the engine is assumed to have the data.
*/
class ManagedVector : public Vector, public DataChangeFirer
{

public:
    /** \brief Construct a ManagedVector from a set of parameters
      * \param[in] params  The description of the ManagedVector
      */
    ManagedVector ( VectorParameters::shared_ptr params );

    /** \brief Construct a view of an AMP vector
      * \param[in] alias  Vector to view
      */
    ManagedVector ( const Vector::shared_ptr alias );

    //! Destructor
    virtual ~ManagedVector ();

    /** \brief  If a vector has multiple views to multiple external packages
      * associated with it, this will return the barest version of the vector
      * \return A vector with the fewest views associated with it.
      * \details  A ManagedVector must have an engine and it may have data.
      * If it has an engine with no data, then the engine has must have data.
      * If the engine can be cast to a ManagedVector, it is and getRootVector
      * is called recursively.
      */
    Vector::shared_ptr  getRootVector();

    /** \brief  Return the engine associated with this ManagedVector
    * \return The engine
    */
    VectorEngine::shared_ptr  getVectorEngine();
    std::string type () const;
    virtual Vector::const_iterator begin() const;
    virtual Vector::const_iterator end() const;
    virtual Vector::iterator begin();
    virtual Vector::iterator end();

    virtual Vector::shared_ptr  subsetVectorForVariable ( const Variable::shared_ptr &name );
    virtual Vector::const_shared_ptr  constSubsetVectorForVariable ( const Variable::shared_ptr &name ) const;
    virtual size_t  numberOfDataBlocks () const;
    virtual size_t  sizeOfDataBlock ( size_t i ) const;
    virtual void copyVector ( const Vector::const_shared_ptr &src_vec  );
    virtual void swapVectors ( Vector &other );
    virtual void aliasVector ( Vector &other );

    virtual bool isAnAliasOf ( Vector &rhs );
    virtual bool isAnAliasOf ( Vector::shared_ptr rhs );
    using Vector::cloneVector;
    virtual boost::shared_ptr<Vector>  cloneVector ( const Variable::shared_ptr name ) const;
    virtual boost::shared_ptr<ParameterBase>  getParameters () ;

    virtual boost::shared_ptr<ManagedVectorParameters>  getManagedVectorParameters () ;

    virtual size_t getLocalSize() const;
    virtual size_t getGlobalSize() const;

    virtual void getValuesByGlobalID ( int numVals , size_t *ndx , double *vals ) const;
    virtual void getLocalValuesByGlobalID ( int numVals , size_t *ndx , double *vals ) const;
    virtual void getGhostValuesByGlobalID ( int numVals , size_t *ndx , double *vals ) const;
    virtual void setValuesByGlobalID(int i, size_t * , const double *val);
    virtual void setLocalValuesByGlobalID(int i, size_t * , const double *val);
    virtual void setGhostValuesByGlobalID(int i, size_t * , const double *val);

    virtual void setToScalar(double alpha);
    virtual void scale(double alpha, const VectorOperations &x);
    virtual void scale(double alpha);
    virtual void add(const VectorOperations &x, const VectorOperations &y);
    virtual void subtract(const VectorOperations &x, const VectorOperations &y);
    virtual void multiply(const VectorOperations &x, const VectorOperations &y);
    virtual void divide(const VectorOperations &x, const VectorOperations &y);
    virtual void reciprocal(const VectorOperations &x);
    virtual void linearSum(double alpha, const VectorOperations &x,
      double beta, const VectorOperations &y);
    virtual void axpy(double alpha, const VectorOperations &x, const VectorOperations &y);
    virtual void axpby(double alpha, double beta, const VectorOperations &x);
    virtual void abs(const VectorOperations &x);
    virtual double min(void) const;
    virtual double max(void) const;
    virtual void setRandomValues(void);
    virtual void setValuesByLocalID(int i, size_t * , const double *val);
    virtual void addValuesByLocalID(int i, size_t * , const double *val);
    virtual void addLocalValuesByGlobalID(int i, size_t * , const double *val);
    virtual void putRawData ( double *in );
    virtual void copyOutRawData ( double **in );
    double L1Norm(void) const;
    double L2Norm(void) const;
    double maxNorm(void) const;
    using Vector::dot;
    double dot(const VectorOperations &x) const;
    virtual UpdateState  getUpdateStatus() const;
    virtual void  setUpdateStatus( UpdateState state );

protected:

    virtual void  selectInto ( const VectorSelector & , shared_ptr );
    virtual void  constSelectInto ( const VectorSelector &criterion , Vector::shared_ptr vector ) const;


    /**\brief  A method that is called whenever data changes.  This fires
         triggers that may have been registered with DataChangeFirer
         */
    virtual void dataChanged ();

    /**\brief  The buffer used to store data
     */
    VectorEngine::BufferPtr    d_vBuffer;
    /**\brief  The engine to act on the buffer
     */
    VectorEngine::shared_ptr     d_Engine;
    /**\brief  The parameters used to create this vector
     */
     boost::shared_ptr<ManagedVectorParameters>  d_pParameters;

    /**\brief  Function that returns a pointer to a managed vector
     */
    virtual ManagedVector *getNewRawPtr () const = 0;
    virtual void *getRawDataBlockAsVoid ( size_t i );
    virtual const void *getRawDataBlockAsVoid ( size_t i ) const;

    virtual void addCommunicationListToParameters ( CommunicationList::shared_ptr comm );

private:

    ManagedVector ();

};


}
}

#include "ManagedVector.inline.h"

#endif
