#ifndef included_AMP_Property
#define included_AMP_Property

#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/MultiVector.h"
#include <algorithm>
#include <limits>
#include <map>
#include <map>
#include <string>
#include <valarray>
#include <vector>

namespace AMP {
namespace Materials {
/**
 * \namespace Materials
 * The materials design envisions a material property as a function.
 * These classes encapsulate the notion of functions that may
 * or may not have parameters that are remembered by the function
 * so you don't have to pass them in at every function call.
 * It also includes fields for sources such as a journal reference,
 * and a name.
 */

/**
 * \class          Property
 * \brief          Provides material properties of scalar type.
 * \tparam        Number  Precision to be used (float, double). Currently, only double is used.
 *
 * A Property class may provide only one of eval(), evalVector() or evalTensor(). It can not
 * provide both scalar and tensor
 */
template <class Number>
class Property {
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
     */
    Property( const std::string &name    = std::string( "NotDefined" ),
              const std::string &source  = std::string( "None" ),
              const Number *params       = NULL,
              const unsigned int nparams = 0,
              const std::string *args    = NULL,
              const unsigned int nargs   = 0,
              const Number ranges[][2]   = NULL )
        : d_name( name ),
          d_source( source ),
          d_params( std::valarray<Number>( params, nparams ) ),
          d_nparams( nparams ),
          d_n_arguments( nargs ),
          d_arguments( nargs ),
          d_defaults( nargs ),
          d_ranges( nargs, std::vector<Number>( 2 ) ),
          d_variableNumberParameters( false )
    {
        if ( args != nullptr ) {
            for ( size_t i = 0; i < d_n_arguments; i++ ) d_arguments[i] = args[i];
        }
        for ( size_t i = 0; i < d_n_arguments; i++ ) {
            d_argToIndexMap.insert( std::pair<std::string, size_t>( d_arguments[i], i ) );
        }
        if ( ranges == NULL && d_n_arguments > 0 ) {
            AMP_INSIST( false, "argument ranges not set" );
        }
        for ( size_t i = 0; i < d_n_arguments; i++ ) {
            for ( size_t j = 0; j < 2; j++ ) {
                d_ranges[i][j] = ranges[i][j];
            }
        }
        for ( size_t i = 0; i < d_n_arguments; i++ ) {
            d_defaults[i] = d_ranges[i][0];
        }
        d_defaultsAreSet = true;
    }

    /**
     * Destructor
     */
    virtual ~Property() {}

    /** return name of property */
    std::string get_name() { return d_name; }

    /** return source reference */
    std::string get_source() { return d_source; }

    /** return property parameters */
    std::valarray<Number> get_parameters() { return d_params; }

    /**
     * \brief		   set the property parameters
     * \param[in]	   params the new parameters
     * \param[in]	   nparams the number of new parameters
     */
    void set_parameters( const Number *params, const unsigned int nparams )
    {
        AMP_INSIST( d_nparams == nparams, "new parameters must be same in number as old" );
        d_params = std::valarray<Number>( params, nparams );
    }

    /**
     * \brief          changing number of parameters allowed
     */
    bool variable_number_parameters() { return d_variableNumberParameters; }

    /**
     * \brief		   set the property parameters
     * \param[in]	   params the new parameters
     * \param[in]	   nparams the number of new parameters
     */
    virtual void set_parameters_and_number( const Number *params, const unsigned int nparams )
    {
        AMP_INSIST( d_variableNumberParameters,
                    "changing number of parameters for this property not allowed" );
        d_params.resize( nparams );
        for ( size_t i = 0; i < nparams; i++ ) d_params[i] = params[i];
        d_nparams                                          = nparams;
    }

    /** return the names of the arguments to eval */
    std::vector<std::string> get_arguments() { return d_arguments; }

    /** return the number of arguments to eval */
    unsigned int get_number_arguments() { return d_n_arguments; }

    /** get the defaults */
    std::vector<Number> get_defaults() { return d_defaults; }

    /** set the defaults */
    void set_defaults( std::vector<Number> defaults )
    {
        AMP_INSIST( defaults.size() == d_n_arguments, "incorrect number of defaults specified" );
        d_defaults = defaults;
    }

    //! get ranges for all arguments used in this material
    virtual std::vector<std::vector<Number>> get_arg_ranges() { return d_ranges; }

    //! determine if a string is an argument
    bool is_argument( const std::string &argname )
    {
        std::map<std::string, size_t>::iterator it = d_argToIndexMap.find( argname );
        if ( it == d_argToIndexMap.end() ) return false;
        return true;
    }

    //! get range for a specific argument
    virtual std::vector<Number> get_arg_range( const std::string &argname )
    {
        std::map<std::string, size_t>::iterator it = d_argToIndexMap.find( argname );
        // AMP_INSIST(it != d_argToIndexMap.end(), std::string("argument ")+argname+std::string("is
        // invalid"));
        if ( it == d_argToIndexMap.end() ) {
            std::vector<Number> infinite_range( 2 );
            infinite_range[0] = -std::numeric_limits<double>::max();
            infinite_range[1] = std::numeric_limits<double>::max();
            return infinite_range;
        }
        size_t index = it->second;
        return d_ranges[index];
    }

    //! determine if a value is within range or not
    bool in_range( const std::string &argname, const Number value );

    //! determine if a set of values are all within range or not
    template <class INPUT_VTYPE>
    bool in_range( const std::string &argname, const INPUT_VTYPE &values );

    //! determine if a set of sets of values are all within range or not
    template <class INPUT_VTYPE>
    bool in_range( const std::map<std::string, AMP::shared_ptr<INPUT_VTYPE>> &values );

    //! set the translation table between property arguments and AMP::Multivector entries
    void set_translator( const std::map<std::string, std::string> &xlator )
    {
        // assure incoming map has correct keys
        for ( std::map<std::string, std::string>::const_iterator iter = xlator.begin();
              iter != xlator.end();
              ++iter ) {
            AMP_ASSERT( std::find( d_arguments.begin(), d_arguments.end(), iter->first ) !=
                        d_arguments.end() );
        }
        d_translator = xlator;
    }

    //! get the translation table between property arguments and AMP::Multivector entries
    std::map<std::string, std::string> get_translator() { return d_translator; }

    //! converts AMP::MultiVector to a map of pointers to AMP::Vectors based on argument names
    std::map<std::string, AMP::shared_ptr<AMP::LinearAlgebra::Vector>>
    make_map( const AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> &args );

    //! indicator for scalar evaluator
    virtual bool isScalar() { return true; }

    //! indicator for vector evaluator
    virtual bool isVector() { return false; }

    //! indicator for tensor evaluator
    virtual bool isTensor() { return false; }

    //! set auxiliary data
    void setAuxiliaryData( const std::string &key, const double val );

    //! set auxiliary data
    void setAuxiliaryData( const std::string &key, const int val );

    //! set auxiliary data
    void setAuxiliaryData( const std::string &key, const std::string &val );

    //! get auxiliary data
    void getAuxiliaryData( const std::string &key, double &val );

    //! get auxiliary data
    void getAuxiliaryData( const std::string &key, int &val );

    //! get auxiliary data
    void getAuxiliaryData( const std::string &key, std::string &val );

protected:
    std::string d_name;             //!< should be unique
    std::string d_source;           //!< journal or report reference: from where did model come?
    std::valarray<Number> d_params; //!< parameters
    unsigned int d_nparams;         //!< number of parameters
    unsigned int d_n_arguments;     //!< number of arguments to the eval function
    std::vector<std::string> d_arguments;          //!< names of the arguments to the eval function
    std::vector<Number> d_defaults;                //!< default values of arguments to eval function
    bool d_defaultsAreSet;                         //!< indicates defaults have been set
    std::vector<std::vector<Number>> d_ranges;     //!< allowed ranges of arguments
    std::map<std::string, size_t> d_argToIndexMap; //!< connects argument names to their indices
    std::map<std::string, std::string> d_translator; //!< standard names to multivector names
    bool d_variableNumberParameters;                 //!< can change number of parameters

    std::map<std::string, double> d_AuxiliaryDataDouble;
    std::map<std::string, int> d_AuxiliaryDataInteger;
    std::map<std::string, std::string> d_AuxiliaryDataString;

    //!//!//!//!//!//!//! Evaluators //!//!//!//!//!//!//!

private:
    /* Loops through input vectors, calling the child eval function, returning scalar */
    template <class INPUT_VTYPE, class RETURN_VTYPE>
    void evalvActual( RETURN_VTYPE &r,
                      const std::map<std::string, AMP::shared_ptr<INPUT_VTYPE>> &args );

public:
    /**
     * scalar evaluation function for a single argument set
     * \param args list of argument values, in correct order, given by  get_arguments()
     * \return scalar value of property
     */
    virtual Number eval( std::vector<Number> &args ) = 0;

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r vector of return values
     *  \param args map of vectors of arguments, indexed by strings which are members of
     * get_arguments()
     *
     *  The  \a args  parameter need not contain all the members of  get_arguments()  as indices.
     *  Arguments left out will have values supplied by the entries in  get_defaults() .
     *  Sizes of  \a r  and \a args["name"] must match. Members of
     *  \a args  indexed by names other than those in  get_arguments()  are ignored.
     */
    void virtual evalv( std::vector<Number> &r,
                        const std::map<std::string, AMP::shared_ptr<std::vector<Number>>> &args );

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r AMP vector of return values
     *  \param args map of AMP vectors of arguments, indexed by strings which are members of
     * get_arguments()
     *
     *  The  \a args  parameter need not contain all the members of  get_arguments()  as indices.
     *  Arguments left out will have values supplied by the entries in  get_defaults() .
     *  Sizes of  \a r  and \a args["name"] must match. Members of
     *  \a args  indexed by names other than those in  get_arguments()  are ignored.
     *  The list {args["name-1"][i], ..., args["name-n"][i]} will be passed to eval() and the k-j-th
     * result
     *  returned in (*r[k][j])[i].
     */
    void virtual evalv(
        AMP::shared_ptr<AMP::LinearAlgebra::Vector> &r,
        const std::map<std::string, AMP::shared_ptr<AMP::LinearAlgebra::Vector>> &args );

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r AMP vector of return values
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
    void virtual evalv( AMP::shared_ptr<AMP::LinearAlgebra::Vector> &r,
                        const AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> &args );
};

template <>
std::map<std::string, AMP::shared_ptr<AMP::LinearAlgebra::Vector>>
Property<double>::make_map( const AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> &args );

template <>
void Property<double>::evalv(
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> &r,
    const std::map<std::string, AMP::shared_ptr<AMP::LinearAlgebra::Vector>> &args );

template <>
void Property<double>::evalv( AMP::shared_ptr<AMP::LinearAlgebra::Vector> &r,
                              const AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> &args );

} // namespace Materials
} // namespace AMP

#include "Property.i.h"

#endif
