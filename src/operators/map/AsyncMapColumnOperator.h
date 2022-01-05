#ifndef included_AMP_AsyncMapColumnOperator
#define included_AMP_AsyncMapColumnOperator

#include "AMP/operators/AsynchronousColumnOperator.h"
#include "AMP/operators/AsynchronousColumnOperatorParameters.h"

#include "AMP/mesh/Mesh.h"

namespace AMP::Operator {

/** \brief   AsyncMapColumnOperator does not require parameters beyond those of
 * AsynchronousColumnOperator
 */
typedef AsynchronousColumnOperatorParameters AsyncMapColumnOperatorParameters;


/** \brief  A column of map operators used for Mesh to Mesh interactions through a mapping.
 * \details  Mappings often require overlapping communication to prevent serial threads
 * from arising during the computation.  Mappings often require results to be stored
 * in auxiliary or "frozen" vectors.  Therefore, the AsyncMapColumnOperator supports a
 * setVector method that will call setVector on all maps in the column.
 *
 * \see AsyncMapOperator
 */
class AsyncMapColumnOperator : public AsynchronousColumnOperator
{
public:
    //! Default constructor
    explicit AsyncMapColumnOperator();

    //! Empty constructor
    explicit AsyncMapColumnOperator( std::shared_ptr<const OperatorParameters> params );

    //! Return the name of the operator
    std::string type() const override { return "AsyncMapColumnOperator"; }

    /** \brief  Call setVector on all vectors in the column
     * \param[in] p  The auxiliary or "frozen" vector to store results in
     */
    void setVector( AMP::LinearAlgebra::Vector::shared_ptr p );

    //!  Returns the frozen vector
    virtual AMP::LinearAlgebra::Vector::shared_ptr getFrozenVector() { return d_OutputVector; }

    void append( std::shared_ptr<Operator> op ) override;

    // Overload the apply operator to include makeConsistent
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /** \brief  A factory method.
     * \return  A column of map operators of type MAP_TYPE
     * \details  In order to use this builder, the MAP_TYPE must have 3 things:  a static member
     *    called CommTagBase, a static method:  bool validMapType ( const std::string & ), and
     *    a typedef for the parameters of the class called Parameters.
     * \tparam  MAP_TYPE  The type of map to create
     * \param[in]  manager     The multi-mesh to search for maps to add to the column
     * \param[in]  database    The global database of the simulation
     *
     */
    template<typename MAP_TYPE>
    static std::shared_ptr<AsyncMapColumnOperator> build( AMP::Mesh::Mesh::shared_ptr manager,
                                                          std::shared_ptr<AMP::Database> database );

    // Function to determine if a makeConsistentSet is required
    virtual bool requiresMakeConsistentSet();

private:
    // Frozen vector for the output results
    AMP::LinearAlgebra::Vector::shared_ptr d_OutputVector;

    // Function to create the databases for the individual maps
    static std::vector<std::shared_ptr<AMP::Database>>
    createDatabases( std::shared_ptr<AMP::Database> database );
};

} // namespace AMP::Operator

#include "AsyncMapColumnOperator.tmpl.h"

#endif
