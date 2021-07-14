#ifndef included_AMP_TensorProperty
#define included_AMP_TensorProperty

#include "Property.h"

namespace AMP {
namespace Materials {

class TensorProperty : public Property
{
public:
    /**
     * Constructor
     * \param name      name of property (required)
     * \param source    literature reference for model and notes
     * \param params    default parameter values
     * \param args      names of arguments
     * \param ranges    ranges of arguments
     * \param dimensions the first and second dimensions of return value tensor
     */
    TensorProperty(
        std::string name,
        std::string source                        = std::string( "None" ),
        std::vector<double> params                = std::vector<double>(),
        std::vector<std::string> args             = std::vector<std::string>(),
        std::vector<std::array<double, 2>> ranges = std::vector<std::array<double, 2>>(),
        std::vector<size_t> dimensions            = { 2, 1 } );

    /**
     * Destructor
     */
    virtual ~TensorProperty() {}

    /// get dimensions of evalv return value tensor
    std::vector<size_t> get_dimensions() { return d_dimensions; }

    /// ok to change dimensions
    bool variable_dimensions() { return d_variableDimensions; }

    /// set dimensions of evalv return value tensor
    void set_dimensions( const std::vector<size_t> &dimensions )
    {
        AMP_INSIST( d_variableDimensions, "can not change dimensions for this property" );
        d_dimensions = dimensions;
        AMP_INSIST( d_dimensions.size() == 2, "there must be two dimensions" );
        AMP_INSIST( d_dimensions[0] > 0, "must have first return tensor dimension > 0" );
        AMP_INSIST( d_dimensions[1] > 0, "must have second return tensor dimension > 0" );
    }

    /// indicator for scalar evaluator
    bool isScalar() override { return false; }

    /// indicator for tensor evaluator
    bool isTensor() override { return true; }

protected:
    std::vector<size_t> d_dimensions; ///< dimensions of return value tensor
    bool d_variableDimensions;        ///< true if ok to change dimensions

    ///////////////////// Evaluators /////////////////////

private:
    /* Loops through input vectors, calling the child eval function, returning tensor results */
    template<class INPUT_VTYPE, class RETURN_VTYPE>
    void evalvActual( std::vector<std::vector<std::shared_ptr<RETURN_VTYPE>>> &r,
                      const std::map<std::string, std::shared_ptr<INPUT_VTYPE>> &args );

public:
    /**
     * tensor evaluation function for a single argument set
     * \param args list of argument values, in correct order, given by  get_arguments()
     * \return tensor of property values with dimensions get_dimensions()
     */
    virtual std::vector<std::vector<double>> evalTensor( std::vector<double> &args ) = 0;

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r tensor vector of return values
     *  \param args map of vectors of arguments, indexed by strings which are members of
     * get_arguments()
     *
     *  The  \a r parameter must have dimensions get_dimensions().
     *  The  \a args  parameter need not contain all the members of  get_arguments()  as indices.
     *  Arguments left out will have values supplied by the entries in  get_defaults() .
     *  Sizes of  *r[i][j]  and \a args["name"] must match. Members of
     *  \a args  indexed by names other than those in  get_arguments()  are ignored.
     *  The list {args["name-1"][i], ..., args["name-n"][i]} will be passed to eval() and the j-th
     * result
     *  returned in (*r[j])[i].
     */
    virtual void evalv( std::vector<std::vector<std::shared_ptr<std::vector<double>>>> &r,
                        const std::map<std::string, std::shared_ptr<std::vector<double>>> &args );

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r tensor of AMP vectors of return values
     *  \param args map of AMP vectors of arguments, indexed by strings which are members of
     * get_arguments()
     *
     *  The  \a args  parameter need not contain all the members of  get_arguments()  as indices.
     *  Arguments left out will have values supplied by the entries in  get_defaults() .
     *  Sizes of  \a *r[i][j]  and \a args["name"] must match. Members of
     *  \a args  indexed by names other than those in  get_arguments()  are ignored.
     *  The list {args["name-1"][i], ..., args["name-n"][i]} will be passed to eval() and the k-j-th
     * result
     *  returned in (*r[k][j])[i].
     */
    virtual void
    evalv( std::vector<std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>>> &r,
           const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args );

    /** Wrapper function that calls evalvActualVector for each argument set
     *  \param r tensor of AMP vectors of return values
     *  \param args AMP multivector of arguments
     *
     *  Before this function is used, a translation table must be assigned by means of the
     * set_translator() function
     *  which gives the correspondence between entries in get_arguments() and the \a args
     * multivector.
     *  Upon invocation, the \a args parameter is converted to a map of AMP vectors via make_map()
     * and passed to another
     * version of evalv.
     */
    virtual void evalv( std::vector<std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>>> &r,
                        const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args );

    // disable scalar evaluator
    double eval( std::vector<double> & ) override
    {
        AMP_ERROR( "cannot use scalar evaluator from tensor property" );
        return 0;
    }

    // disable scalar evaluator
    void evalv( std::vector<double> &,
                const std::map<std::string, std::shared_ptr<std::vector<double>>> & ) override
    {
        AMP_ERROR( "cannot use scalar evaluator from tensor property" );
    }

    // disable scalar evaluator
    void
    evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &,
           const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> & ) override
    {
        AMP_ERROR( "cannot use scalar evaluator from tensor property" );
    }

    // disable scalar evaluator
    void evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &,
                const std::shared_ptr<AMP::LinearAlgebra::MultiVector> & ) override
    {
        AMP_ERROR( "cannot use scalar evaluator from tensor property" );
    }
};


} // namespace Materials
} // namespace AMP


#endif /* included_AMP_TensorProperty */
