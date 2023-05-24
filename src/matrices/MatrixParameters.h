#ifndef included_AMP_MatrixParameters
#define included_AMP_MatrixParameters

#include "AMP/vectors/Vector.h"


namespace AMP::Discretization {
class DOFManager;
}


namespace AMP::LinearAlgebra {


/** \class MatrixParameters
 * \brief  A class used to hold basic parameters for a matrix
 */
class MatrixParameters
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
    explicit MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> left,
                               std::shared_ptr<AMP::Discretization::DOFManager> right,
                               const AMP_MPI &comm );

    //! Deconstructor
    virtual ~MatrixParameters(){};

    //! Return the local number of rows
    size_t getLocalNumberOfRows() const;

    //! Return the local number of columns
    size_t getLocalNumberOfColumns() const;

    //! Return the global number of rows
    size_t getGlobalNumberOfRows() const;

    //! Return the global number of columns
    size_t getGlobalNumberOfColumns() const;

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

    /** \brief Return a vector of all columns with non-zero entries for this process
     * \return a vector of global column ids
     */
    std::vector<size_t> &getColumns() { return d_vColumns; }

    /** \brief  Add columns to a description
     * \param[in] cols  The column ids
     */
    void addColumns( const std::vector<size_t> &cols );

    //!  Get the DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a
    //!  left vector )
    std::shared_ptr<AMP::Discretization::DOFManager> getLeftDOFManager();

    //!  Get the DOFManager for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a
    //!  right vector )
    std::shared_ptr<AMP::Discretization::DOFManager> getRightDOFManager();

    //!  Get the communicator for the matrix
    AMP::AMP_MPI &getComm();

    //!  The communication list of a left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a
    //!  left vector )
    std::shared_ptr<CommunicationList> d_CommListLeft;

    //!  The communication list of a right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a
    //!  right vector )
    std::shared_ptr<CommunicationList> d_CommListRight;

    //!  The variable for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left
    //!  vector )
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_VariableLeft;

    //!  The variable for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right
    //!  vector )
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_VariableRight;

protected:
    MatrixParameters(){};

    // The DOFManager for the left vector ( may be null )
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFManagerLeft;

    // The DOFManager for the right vector ( may be null )
    std::shared_ptr<AMP::Discretization::DOFManager> d_DOFManagerRight;

    //!  The number of nonzeros per row of the matrix
    std::vector<int> d_vEntriesPerRow;

    //!  The set of columns with non-zero entries this processor has
    std::vector<size_t> d_vColumns;

    // The comm of the matrix
    AMP_MPI d_comm;
};
} // namespace AMP::LinearAlgebra

#endif
