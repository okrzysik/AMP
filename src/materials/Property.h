#ifndef included_AMP_Property
#define included_AMP_Property

#include "AMP/utils/Units.h"
#include "AMP/utils/UtilityMacros.h"

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>


// Foward declare classes
namespace AMP::LinearAlgebra {
class Vector;
class MultiVector;
} // namespace AMP::LinearAlgebra


namespace AMP::Materials {


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
 *
 * A Property class may provide only one of eval(), evalVector() or evalTensor(). It can not
 * provide both scalar and tensor
 */
class Property
{
public:
    /**
     * Constructor
     * \param name      name of property (required)
     * \param source    literature reference for model and notes
     * \param params    default parameter values
     * \param args      names of arguments
     * \param units     Optional units for each argument
     */
    Property( std::string name,
              std::string source                        = "None",
              std::vector<double> params                = {},
              std::vector<std::string> args             = {},
              std::vector<std::array<double, 2>> ranges = {},
              std::vector<Units> units                  = {} );

    /**
     * Destructor
     */
    virtual ~Property() {}

    /** return name of property */
    inline std::string get_name() const { return d_name; }

    /** return source reference */
    inline std::string get_source() const { return d_source; }

    /** return property parameters */
    inline const std::vector<double> &get_parameters() const { return d_params; }

    /**
     * \brief		   set the property parameters
     * \param[in]	   params the new parameters
     * \param[in]	   nparams the number of new parameters
     */
    inline void set_parameters( std::vector<double> params )
    {
        AMP_INSIST( d_params.size() == params.size(),
                    "new parameters must be same in number as old" );
        d_params = std::move( params );
    }

    /**
     * \brief          changing number of parameters allowed
     */
    inline bool variable_number_parameters() { return d_variableNumberParameters; }

    /**
     * \brief		   set the property parameters
     * \param[in]	   params the new parameters
     * \param[in]	   nparams the number of new parameters
     */
    virtual void set_parameters_and_number( std::vector<double> params );

    /** return the names of the arguments to eval */
    inline const std::vector<std::string> &get_arguments() const { return d_arguments; }

    /** return the number of arguments to eval */
    inline unsigned int get_number_arguments() { return d_arguments.size(); }

    /** get the defaults */
    inline const std::vector<double> &get_defaults() { return d_defaults; }

    /** set the defaults */
    inline void set_defaults( std::vector<double> defaults )
    {
        AMP_INSIST( defaults.size() == d_arguments.size(),
                    "incorrect number of defaults specified" );
        d_defaults = defaults;
    }

    //! determine if a string is an argument
    bool is_argument( const std::string &argname );

    //! indicator for scalar evaluator
    virtual bool isScalar() { return true; }

    //! indicator for vector evaluator
    virtual bool isVector() { return false; }

    //! indicator for tensor evaluator
    virtual bool isTensor() { return false; }

    // converts AMP::MultiVector to a map of pointers to AMP::Vectors based on argument names
    std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>>
    make_map( const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
              const std::map<std::string, std::string> &translator ) const;


public: // Functions dealing with the ranges of the arguments
    //! get units for all arguments used in this material
    std::vector<Units> get_arg_units() { return d_argUnits; }

    //! get ranges for all arguments used in this material
    std::vector<std::array<double, 2>> get_arg_ranges() { return d_ranges; }

    //! get range for a specific argument
    std::array<double, 2> get_arg_range( const std::string &argname );

    //! determine if a value is within range or not
    inline bool in_range( const std::string &argname,
                          double value,
                          Units unit      = Units(),
                          bool throwError = false );

    //! determine if a set of values are all within range or not
    template<class INPUT_VTYPE>
    inline bool in_range( const std::string &argname,
                          const INPUT_VTYPE &values,
                          Units unit      = Units(),
                          bool throwError = false );


public: // Functions dealing with auxilliary data
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


public: // Evaluators
    /**
     * scalar evaluation function for a single argument set
     * \param args list of argument values, in correct order, in the correct units, given by
     *    get_arguments()
     * \return scalar value of property
     */
    virtual double eval( const std::vector<double> &args ) = 0;

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r vector of return values
     *  \param args map of vectors of arguments, indexed by strings which are members of
     *     get_arguments()
     *
     *  The  \a args  parameter need not contain all the members of  get_arguments()  as
     * indices. Arguments left out will have values supplied by the entries in  get_defaults() .
     *  Sizes of  \a r  and \a args["name"] must match. Members of
     *  \a args  indexed by names other than those in  get_arguments()  are ignored.
     */
    void evalv( std::vector<double> &r,
                const std::map<std::string, std::shared_ptr<std::vector<double>>> &args );

    /** Wrapper function that calls evalvActual for each argument set
     *  \param r AMP vector of return values
     *  \param args map of AMP vectors of arguments, indexed by strings which are members of
     * get_arguments()
     *
     *  The  \a args  parameter need not contain all the members of  get_arguments()  as
     * indices. Arguments left out will have values supplied by the entries in  get_defaults() .
     *  Sizes of  \a r  and \a args["name"] must match. Members of
     *  \a args  indexed by names other than those in  get_arguments()  are ignored.
     *  The list {args["name-1"][i], ..., args["name-n"][i]} will be passed to eval() and the
     * k-j-th result returned in (*r[k][j])[i].
     */
    void evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
                const std::map<std::string, std::shared_ptr<AMP::LinearAlgebra::Vector>> &args );

    /** Wrapper function that calls evalvActual for each argument set
     *  Upon invocation, the \a args parameter is converted to a map of AMP vectors via
     *     make_map() and passed to another version of evalv.
     *  \param r            AMP vector of return values
     *  \param args         AMP multivector of arguments
     *  \param translator   Optional translator between property arguments and AMP::Multivector
     * entries
     */
    void evalv( std::shared_ptr<AMP::LinearAlgebra::Vector> &r,
                const std::shared_ptr<AMP::LinearAlgebra::MultiVector> &args,
                const std::map<std::string, std::string> &translator = {} );


protected:
    Units d_units;                        //!< default units to return
    std::string d_name;                   //!< should be unique
    std::string d_source;                 //!< reference for source data
    std::vector<double> d_params;         //!< parameters
    std::vector<std::string> d_arguments; //!< names of the arguments to the eval function
    std::vector<Units> d_argUnits;        //!< default units for the arguments to the eval function
    std::vector<double> d_defaults;       //!< default values of arguments to eval function
    bool d_defaultsAreSet;                //!< indicates defaults have been set
    std::vector<std::array<double, 2>> d_ranges;   //!< allowed ranges of arguments
    std::map<std::string, size_t> d_argToIndexMap; //!< connects argument names to their indices
    bool d_variableNumberParameters;               //!< can change number of parameters

    std::map<std::string, double> d_AuxiliaryDataDouble;
    std::map<std::string, int> d_AuxiliaryDataInteger;
    std::map<std::string, std::string> d_AuxiliaryDataString;


protected:
    // Get the index for the desired argument
    inline int get_arg_index( const std::string &name )
    {
        auto it = d_argToIndexMap.find( name );
        if ( it == d_argToIndexMap.end() )
            return -1;
        return it->second;
    }

    // Loops through input vectors, calling the child eval function, returning scalar
    template<class INPUT_VTYPE, class RETURN_VTYPE>
    void evalvActual( RETURN_VTYPE &r,
                      const std::map<std::string, std::shared_ptr<INPUT_VTYPE>> &args );
};


} // namespace AMP::Materials

#include "Property.i.h"

#endif
