#ifndef included_AMP_ManagedMatrixParameters
#define included_AMP_ManagedMatrixParameters


#include "AMP/matrices/MatrixParameters.h"

#include <set>


namespace AMP {
namespace LinearAlgebra {


/** \class     ManagedMatrixParameters
 * \brief  A class used to create an Epetra matrix
 */
class ManagedMatrixParameters : public MatrixParameters
{
public:
    /** \brief Constructor
     * \param[in] left     The DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$y\f$ is a left
     * vector )
     * \param[in] right    The DOFManager for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$x\f$ is a right
     * vector )
     * \param[in] comm     Communicator for the matrix
     */
    explicit ManagedMatrixParameters( AMP::Discretization::DOFManager::shared_ptr left,
                                      AMP::Discretization::DOFManager::shared_ptr right,
                                      const AMP_MPI &comm );


    //! Deconstructor
    virtual ~ManagedMatrixParameters(){};

    /** \brief Return the number of entries in each row
     * \return  An integer array of the number of entries in each
     * local row
     */
    const int *entryList() const;

    /** \brief Return the number of entries in each row
     * \return  An integer array of the number of entries in each
     * local row
     */
    int *entryList();

    /** \brief Set the number of non-zeros in a particular row
     * \param[in] row  The row number
     * \param[in] entries  The number of non-zero entries
     */
    void setEntriesInRow( int row, int entries );

    /** \brief Return the number of non-zero entries in a local row
     * \param[in] i The local row id
     * \return  The number of entries in the row
     */
    int &entriesInRow( int i );

    /** \brief Return the number of non-zero entries in a local row
     * \param[in] i The local row id
     * \return  The number of entries in the row
     */
    int entriesInRow( int i ) const;

    /** \brief Return the maximum number of non-zero entries in a
     * local row.  This is not a global operation.
     * \return  The maximum number of entries in any row on this core.
     */
    int maxEntitiesInRow() const;


    /** \brief  Is the matrix described square
     * \return True if the number of cols = number of rows
     */
    bool isSquare();

    /** \brief  Add columns to a description
     * \param[in] i  The number of columns
     * \param[in] cols  The column ids
     */
    void addColumns( int i, int *cols );

    /** \brief  Add columns to a description
     * \param[in] cols  The column ids
     */
    void addColumns( const std::set<size_t> &cols );

protected:
    //! Constructor -- unimplemented
    ManagedMatrixParameters();

    //! Constructor -- unimplemented
    ManagedMatrixParameters( const ManagedMatrixParameters & );

    //!  The number of nonzeros per row of the matrix
    std::vector<int> d_vEntriesPerRow;

    //!  The set of columns this processor has
    std::set<int> d_sColumns;
};
} // namespace LinearAlgebra
} // namespace AMP

#endif
