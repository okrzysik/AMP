#ifndef included_AMP_ManagedEpetraMatrix
#define included_AMP_ManagedEpetraMatrix

#include <set>

#include "EpetraMatrix.h"
#include "ManagedMatrix.h"
#include "vectors/ManagedDataMap.h"
#include "vectors/trilinos/EpetraVector.h"

#include <Epetra_FECrsMatrix.h>

namespace AMP {
namespace LinearAlgebra {

  /** \class ManagedEpetraMatrixParameters
    * \brief  A class used to create an Epetra matrix
    */
  class ManagedEpetraMatrixParameters : public MatrixParameters ,
                                        public ManagedDataMap
  {
    private:
      boost::shared_ptr < Epetra_Map >   d_eRowMap;
      boost::shared_ptr < Epetra_Map >   d_eColMap;
      AMP_MPI                            d_comm;

    protected:
      /** \brief Constructor -- unimplemented
        */
      ManagedEpetraMatrixParameters ();

      /** \brief Constructor -- unimplemented
        */
      ManagedEpetraMatrixParameters ( const ManagedEpetraMatrixParameters & );

      /** \brief  The number of nonzeros per row of the matrix
        */
      std::vector<int>             d_vEntriesPerRow;

      /** \brief  The set of columns this processor has
        */
      std::set<int>                d_sColumns;

      /** \brief  Total number of columns
        */
      int    d_ColGlobal;

      /** \brief  First column on this core
        */
      int    d_ColBase;

      /** \brief  First row on this core
        */
      int    d_RowBase;

    public:
      /** \brief Convenience typedef
        */
      typedef  ManagedDataMap::iterator              iterator;

      /** \brief Convenience typedef
        */
      typedef  ManagedDataMap::const_iterator        const_iterator;

      /** \brief The communication list of a left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left vector )
        */
      CommunicationList::shared_ptr   d_CommListLeft;
      /** \brief The communication list of a right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
        */
      CommunicationList::shared_ptr   d_CommListRight;

      /** \brief Constructor
        * \param[in] local_size  Number of rows on this core
        * \param[in] global_size Number of total rows 
        * \param[in] first_dof  Global ID of first row
        * \param[in] c Communicator for the matrix
        */
      ManagedEpetraMatrixParameters ( int local_size , int global_size , int first_dof , AMP_MPI c );

      /** \brief Constructor
        * \param[in] row_local_size  Number of rows on this core
        * \param[in] row_global_size Number of rows in the matrix
        * \param[in] row_first_dof   ID of the first row on this core
        * \param[in] col_global_size Number of columns in the matrix
        * \param[in] col_first_dof   ID of the first column on this core
        * \param[in] c The communicator for the matrix
        */
      ManagedEpetraMatrixParameters ( int row_local_size , int row_global_size , int row_first_dof , int col_global_size , int col_first_dof , AMP_MPI c );

      /** \brief Return the number of entries in each row
        * \return  An integer array of the number of entries in each
        * local row
        */
      const int * entryList () const;

      /** \brief Return the number of entries in each row
        * \return  An integer array of the number of entries in each
        * local row
        */
            int * entryList(); 

      /** \brief Set the number of non-zeros in a particular row
        * \param[in] row  The row number
        * \param[in] entries  The number of non-zero entries
        */
      void   setEntriesInRow ( int row , int entries );

      /** \brief Return the number of non-zero entries in a local row
        * \param[in] i The local row id
        * \return  The number of entries in the row
        */
      int  & entriesInRow ( int i );

      /** \brief Return the number of non-zero entries in a local row
        * \param[in] i The local row id
        * \return  The number of entries in the row
        */
      int    entriesInRow ( int i ) const;

      /** \brief Return the maximum number of non-zero entries in a 
        * local row.  This is not a global operation.
        * \return  The maximum number of entries in any row on this core.
        */
      int    maxEntitiesInRow () const;


      /** \brief  Is the matrix described square
        * \return True if the number of cols = number of rows
        */
      bool  isSquare ();

      /** \brief  Get the Epetra_Map for the rows
        * \return  The Epetra_Map
        */
      Epetra_Map      &getEpetraRowMap ();

      /** \brief  Get the Epetra_Map for the columns
        * \return  The Epetra_Map
        */
      Epetra_Map      *getEpetraColMap ();

      /** \brief  Get the Epetra_Map for the rows as a shared pointer
        * \return  The Epetra_Map
        */
      boost::shared_ptr < Epetra_Map >   getEpetraRowMapPtr ();

      /** \brief  Get the Epetra_Map for the columns as a shared pointer
        * \return  The Epetra_Map
        */
      boost::shared_ptr < Epetra_Map >   getEpetraColMapPtr ();

      /** \brief  Get the AMP_MPI comm associated with this description
        * \return  The AMP_MPI object
        */
      AMP_MPI   getEpetraComm ();

      /** \brief  Add columns to a description
        * \param[in] i  The number of columns
        * \param[in] cols  The column ids
        */
      void             addColumns ( int i , int *cols );
  };



  /** \class ManagedEpetraMatrix
    * \brief  A class that wraps an Epetra_CrsMatrix
    * \details  This class stores an Epetra_FECrsMatrix and provides
    * the AMP interface to this matrix.
    */
  class ManagedEpetraMatrix : public EpetraMatrix ,
                              public ManagedMatrix
  {
    protected:
      /** \brief  Parameters used to construct the matrix
        */
      ParametersPtr         d_pParameters;
      
      /** \brief  Unimplemented constructor
        */
      ManagedEpetraMatrix();

      /** \brief  Unimplemented constructor
        */
      ManagedEpetraMatrix ( const ManagedEpetraMatrix &rhs );

      /** \brief  \f$A_{i,j}\f$ storage of off-core data
        */
      std::map<int,std::map<int,double> >  d_OtherData;

      /** \brief  Update data off-core
        */
      void  setOtherData ();

      virtual void multiply ( shared_ptr other_op , shared_ptr &result );

    public:
      /** \brief Constructor
        * \param[in] p  The description of the matrix
        */
      ManagedEpetraMatrix( ParametersPtr p );

      /** \brief Constructor from Epetra_CrsMatrix
        * \param[in]  m  Matrix to wrap
        * \param[in]  dele  If true, this class deletes the matrix
        */
      ManagedEpetraMatrix ( Epetra_CrsMatrix *m , bool dele = false );

      /** \brief Destructor
        */
      virtual ~ManagedEpetraMatrix() {}

      virtual void  createValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values );



      virtual void mult(const Vector::shared_ptr & in, Vector::shared_ptr & out);
      virtual void multTranspose (const Vector::shared_ptr & in, Vector::shared_ptr & out);
      virtual Vector::shared_ptr  extractDiagonal ( Vector::shared_ptr buf = Vector::shared_ptr() );
      virtual void  scale ( double alpha );
      virtual void  axpy ( double alpha , const Matrix &rhs );
      virtual size_t  numRows () { return d_epetraMatrix->NumGlobalRows(); }
      virtual size_t  numColumns () { return d_epetraMatrix->NumGlobalCols(); }
      virtual void  addValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values );
      virtual void  setValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values );

      virtual void  getRowByGlobalID ( int row ,
                                       std::vector<unsigned int> &cols,
                                       std::vector<double>       &values ) const;

      virtual void setScalar ( double );
      virtual void setDiagonal ( const Vector::shared_ptr &in );

      virtual void makeConsistent ();
      virtual double  L1Norm() const;
      virtual Matrix::shared_ptr cloneMatrix () const;
      virtual Vector::shared_ptr  getRightVector ();
      virtual Vector::shared_ptr  getLeftVector ();
      virtual void fillComplete();
  };

}
}

#include "ManagedEpetraMatrix.inline.h"
#endif


