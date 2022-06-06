#ifndef included_AMP_TensorProperty
#define included_AMP_TensorProperty

#include "AMP/materials/Property.h"
#include "AMP/utils/Array.h"


namespace AMP::Materials {

class TensorProperty : public Property
{
public:
    /**
     * Constructor
     * \param name      name of property (required)
     * \param source    literature reference for model and notes
     * \param args      names of arguments
     * \param ranges    ranges of arguments
     * \param dimensions the first and second dimensions of return value tensor
     */
    TensorProperty(
        std::string name,
        std::string source                        = std::string( "None" ),
        std::vector<std::string> args             = std::vector<std::string>(),
        std::vector<std::array<double, 2>> ranges = std::vector<std::array<double, 2>>(),
        const AMP::ArraySize &dimensions          = { 2, 1 } );

    /**
     * Destructor
     */
    virtual ~TensorProperty() {}

    /// get dimensions of evalv return value tensor
    ArraySize get_dimensions() const { return d_dimensions; }

    /// indicator for scalar evaluator
    bool isScalar() const override { return false; }

    /// indicator for tensor evaluator
    bool isTensor() const override { return true; }

protected:
    ArraySize d_dimensions; ///< dimensions of return value tensor

public:
    /**
     * tensor evaluation function for a single argument set
     * \param args list of argument values, in correct order, given by  get_arguments()
     * \return tensor of property values with dimensions get_dimensions()
     */
    virtual AMP::Array<double> evalTensor( const std::vector<double> &args ) const = 0;

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
    void evalv( AMP::Array<std::shared_ptr<std::vector<double>>> &r,
                const std::map<std::string, std::shared_ptr<std::vector<double>>> &args ) const;

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
    void
    evalv( AMP::Array<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
           const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args ) const;

    /** Wrapper function that calls evalvActualVector for each argument set
     *  Upon invocation, the \a args parameter is converted to a map of AMP vectors via make_map()
     *     and passed to another version of evalv.
     *  \param r tensor of AMP vectors of return values
     *  \param args AMP multivector of arguments
     *  \param translator   Optional translator between property arguments and AMP::Multivector
     * entries
     */
    void evalv( AMP::Array<std::shared_ptr<AMP::LinearAlgebra::Vector>> &r,
                const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                const std::map<std::string, std::string> &translator = {} ) const;


public: // Advanced interfaces
    /* Loops through input vectors, calling the child eval function, returning tensor results */
    template<class OUT, class IN = OUT>
    void evalv( AMP::Array<std::shared_ptr<OUT>> &r,
                const std::vector<argumentDataStruct<IN>> &args ) const;

    // disable scalar evaluator
    double eval( const std::vector<double> & ) const override
    {
        AMP_ERROR( "cannot use scalar evaluator from tensor property" );
        return 0;
    }
};


} // namespace AMP::Materials


#endif /* included_AMP_TensorProperty */
