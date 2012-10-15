#ifndef included_AMP_ManagedEpetraMatrixParameters
#define included_AMP_ManagedEpetraMatrixParameters

#include <set>

#include "matrices/MatrixParameters.h"


#include <Epetra_FECrsMatrix.h>

namespace AMP {
namespace LinearAlgebra {


/** \class ManagedEpetraMatrixParameters
  * \brief  A class used to create an Epetra matrix
  */
class ManagedEpetraMatrixParameters : public MatrixParameters
{
public:

    /** \brief Constructor
      * \param[in] left     The DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left vector )
      * \param[in] right    The DOFManager for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
      * \param[in] comm     Communicator for the matrix
      */
    ManagedEpetraMatrixParameters( AMP::Discretization::DOFManager::shared_ptr left, AMP::Discretization::DOFManager::shared_ptr right, AMP_MPI comm );


    //! Deconstructor
    virtual ~ManagedEpetraMatrixParameters() {};

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

private:
    boost::shared_ptr < Epetra_Map >   d_eRowMap;
    boost::shared_ptr < Epetra_Map >   d_eColMap;

protected:
    //! Constructor -- unimplemented
    ManagedEpetraMatrixParameters ();

    //! Constructor -- unimplemented
    ManagedEpetraMatrixParameters ( const ManagedEpetraMatrixParameters & );

    //!  The number of nonzeros per row of the matrix
    std::vector<int>             d_vEntriesPerRow;

    //!  The set of columns this processor has
    std::set<int>                d_sColumns;

};



}
}

#endif


