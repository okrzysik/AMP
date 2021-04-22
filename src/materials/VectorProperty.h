#ifndef included_AMP_VectorProperty
#define included_AMP_VectorProperty

#include "Property.h"

namespace AMP {
namespace Materials {

class VectorProperty : public Property
{
public:
    /**
     * Constructor
     * \param name name of property
     * \param source literature reference for model and notes
     * \param params default parameter values
     * \param nparams number of parameter values
     * \param args names of arguments
     * \param nargs number of arguments
     * \param ranges ranges of arguments
     * \param dimension dimension of return value vector
     */
    VectorProperty( const std::string &name    = std::string( "NotDefined" ),
                    const std::string &source  = std::string( "None" ),
                    const double *params       = nullptr,
                    const unsigned int nparams = 0,
                    const std::string *args    = nullptr,
                    const unsigned int nargs   = 0,
                    const double ranges[][2]   = nullptr,
                    const size_t dimension     = 1 )
        : Property( name, source, params, nparams, args, nargs, ranges ),
          d_dimension( dimension ),
          d_variableDimension( false )
    {
        AMP_INSIST( d_dimension > 0, "must return at least one value" );
    }

    /**
     * Destructor
     */
    virtual ~VectorProperty() {}

    /// return dimension of evalv return value vector
    size_t get_dimension() { return d_dimension; }

    /// ok to change dimension
    bool variable_dimension() { return d_variableDimension; }

    /// supply dimension of evalv return value vector
    void set_dimension( size_t dimension )
    {
        AMP_INSIST( d_variableDimension, "can not change dimension for this property" );
        d_dimension = dimension;
        AMP_INSIST( d_dimension > 0, "must return at least one value" );
    }

    /// indicator for scalar evaluator
    bool isScalar() override { return false; }

    /// indicator for vector evaluator
    bool isVector() override { return true; }

protected:
    size_t d_dimension;       ///< number of return values
    bool d_variableDimension; ///< true if ok to change dimension

    ///////////////////// Evaluators /////////////////////

private:
    /* Loops through input vectors, calling the child eval function, returning vector results */
    template<class INPUT_VTYPE, class RETURN_VTYPE>
    void evalvActual( std::vector<std::shared_ptr<RETURN_VTYPE>> &r,
                      const std::map<std::string, std::shared_ptr<INPUT_VTYPE>> &args );

public:
    /**
     * vector evaluation function for a single argument set
     * \param args list of argument values, in correct order, given by  get_arguments()
     * \return vector of property values of dimension get_dimensions()[0]
     */
    virtual std::vector<double> evalVector( std::vector<double> &args ) = 0;

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r vector of vectors of return values
     *  \param args map of vectors of arguments, indexed by strings which are members of
     * get_arguments()
     *
     *  The  \a r parameter must have size get_dimensions()[0].
     *  The  \a args  parameter need not contain all the members of  get_arguments()  as indices.
     *  Arguments left out will have values supplied by the entries in  get_defaults() .
     *  Sizes of  *r[i]  and \a args["name"] must match. Members of
     *  \a args  indexed by names other than those in  get_arguments()  are ignored.
     *  The list {args["name-1"][i], ..., args["name-n"][i]} will be passed to eval() and the result
     *  returned in r[i].
     */
    virtual void evalv( std::vector<std::shared_ptr<std::vector<double>>> &r,
                        const std::map<std::string, std::shared_ptr<std::vector<double>>> &args );

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r vector of AMP vectors of return values
     *  \param args map of AMP vectors of arguments, indexed by strings which are members of
     * get_arguments()
     *
     *  The  \a args  parameter need not contain all the members of  get_arguments()  as indices.
     *  Arguments left out will have values supplied by the entries in  get_defaults() .
     *  Sizes of  \a *r[i]  and \a args["name"] must match. Members of
     *  \a args  indexed by names other than those in  get_arguments()  are ignored.
     *  The list {args["name-1"][i], ..., args["name-n"][i]} will be passed to eval() and the j-th
     * result
     *  returned in (*r[j])[i].
     */
    virtual void
    evalv( std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
           const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args );

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r vector of AMP vectors of return values
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
    virtual void evalv( std::vector<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                        const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args );

    // disable scalar evaluator
    double eval( std::vector<double> & ) override
    {
        AMP_ERROR( "cannot use scalar evaluator from vector property" );
        return 0;
    }

    // disable scalar evaluator
    void evalv( std::vector<double> &,
                const std::map<std::string, std::shared_ptr<std::vector<double>>> & ) override
    {
        AMP_ERROR( "cannot use scalar evaluator from vector property" );
    }

    // disable scalar evaluator
    void
    evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &,
           const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> & ) override
    {
        AMP_ERROR( "cannot use scalar evaluator from vector property" );
    }

    // disable scalar evaluator
    void evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &,
                const std::shared_ptr<AMP::LinearAlgebra::MultiVector> & ) override
    {
        AMP_ERROR( "cannot use scalar evaluator from vector property" );
    }
};


} // namespace Materials
} // namespace AMP


#endif /* included_AMP_VectorProperty */
