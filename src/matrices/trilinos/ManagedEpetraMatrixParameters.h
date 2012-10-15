#ifndef included_AMP_ManagedEpetraMatrixParameters
#define included_AMP_ManagedEpetraMatrixParameters

#include <set>

#include "matrices/ManagedMatrix.h"
#include "matrices/trilinos/EpetraMatrix.h"
#include "vectors/trilinos/EpetraVector.h"
#include "discretization/DOF_Manager.h"

#include <Epetra_FECrsMatrix.h>

namespace AMP {
namespace LinearAlgebra {


/** \class ManagedEpetraMatrixParameters
  * \brief  A class used to create an Epetra matrix
  */
class ManagedEpetraMatrixParameters : public MatrixParameters
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

      //!  The communication list of a left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left vector )
      CommunicationList::shared_ptr   d_CommListLeft;
      //!  The communication list of a right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
      CommunicationList::shared_ptr   d_CommListRight;
      //!  The DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left vector )
      AMP::Discretization::DOFManager::shared_ptr   d_DOFManagerLeft;
      //!  The DOFManager for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
      AMP::Discretization::DOFManager::shared_ptr   d_DOFManagerRight;
      //!  The variable for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left vector )
      AMP::LinearAlgebra::Variable::shared_ptr   d_VariableLeft;
      //!  The variable for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
      AMP::LinearAlgebra::Variable::shared_ptr   d_VariableRight;

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
      Epetra_Map      &getEpetraColMap ();

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



}
}

#endif


