#ifndef included_AMP_RNG_h
#define included_AMP_RNG_h

#include "utils/shared_ptr.h"
#include <cstdlib>

namespace AMP {


//! Parameters used to instantiate a pseudorandom number generator
class RNGParameters
{
public:
    /**\brief Flag to let the RNG know if you want to provide a seed or use a global seed
     */
    enum class RNGOptions { USE_GLOBAL_SEED, USE_PARAMS_SEED };

    /**\brief Shorthand for shared pointer to RNGParameters
     */
    typedef AMP::shared_ptr<RNGParameters> shared_ptr;

    /**\brief Seed to use when creating an RNG
     */
    size_t d_Seed;

    /**\brief Rank of the RNG to use
     *\details  It is common for simulations to require multiple i.i.d. streams.
     * This parameter will allow the user to select which RNG to create
     */
    size_t d_Rank;

    /**\brief Which seed should be used
     */
    RNGOptions d_WhichSeed;

    /**\brief Constructor.
     */
    RNGParameters( RNGOptions o = RNGOptions::USE_GLOBAL_SEED, size_t rank = 0, size_t seed = 0 );
};


/**\brief A pseudorandom number stream
 *\details This class implements a parallel pseudorandom number stream.  Given a seed
 * and rank, this class will provide a stream of pseudorandom bits, integers, or doubles.
 *
 * USE OF THIS CLASS IS DANGEROUS FOR SIMULATION.  THIS WRAPS THE C FUNCTION RAND() AND
 * DOES NOT GENERATE STREAMS OF SUFFICIENT INDEPENDENCE FOR USE WITH SIMULATION.  USE
 * A DERIVED CLASS FOR BEST RESULTS.
 */
class RNG
{
protected:
    /**\brief  Constant used in the generation of random doubles
     */
    static double d_SizeTDivisor;

    /**\brief  A global seed used for convenience
     */
    static size_t d_Seed;

    /**\brief  Parameters used to construct this class
     */
    RNGParameters::shared_ptr d_Params;

public:
    /**\brief Shorthand for shared pointer to RNG
     */
    typedef AMP::shared_ptr<RNG> shared_ptr;

    /**\brief Initialization function to be called at program start
     *\details  Computes the static constants
     */
    static void initialize( size_t seed = 0 );

    /**\brief Constructor
     * \param[in] params  Description of parameters to create RNG class
     *\details  This calls srand() with the chosen seed plus the rank.  THIS IS A COMPLETELY
     * INADEQUATE RNG.
     */
    explicit RNG( RNGParameters::shared_ptr params );

    /**\brief Destructor
     */
    virtual ~RNG() {}

    /**\brief Fill a buffer with random bits
       \param[in] buf   Buffer to fill
       \param[in] len   Size of buffer in bytes
     */
    virtual void fillBuffer( void *buf, size_t len );

    /**\brief Return the next integer in (low, low+1 ,...,high-1)
     *\param[in] low The smallest integer that may be returned
     *\param[in] high  One more than the largest integer that may be returned
     */
    virtual int nextInt( int low, int high );

    /**\brief Return the next double in [low,high) using a simple algorithm
     *\param[in] low  The smallest double that may be returned
     *\param[in] high  The supremeum of numbers that may be returned.
     *\details  The value is determined by dividing a random size_t with
     * one more than the largest size_t.  d_SizeTDivisor is determined by
     * approximating \f$\epsilon_{\mathit{mach}}\f$ and using this to
     * compute the next largest representable double.
     */
    virtual double nextDouble( double low, double high );

    /**\brief Return a new RNG with a different rank.
     *\param[in]  new_rank  New rank of cloned RNG.
     */
    virtual RNG::shared_ptr cloneRNG( size_t new_rank );
};


/**\brief  A convenience class for using pseudorandom streams
  *\tparam  T  The type to cast as for random values
  *\details  This class can be used with any P.O.D. type to generate a random variable
  * that can be used as a regular variable.  For instance, \code
  void performMath ( RNG::shared_ptr stream )
  {
    RandomVariable<int>     randomInt ( 0 , 100 , stream );
    RandomVariable<double>  randomDouble ( 0. , 1. , stream );
    double  array[100];

    // Assign random values to random spots in an array
    for ( size_t i = 0 ; i != 50 ; i++ )
      array[randomInt] = randomDouble;
  }
  \endcode
  */
template<typename T>
class RandomVariable
{
public:
    /**\brief Typedef of the type of the random variable
     */
    typedef T type;

private:
    type d_Low;
    type d_Size;
    RNG::shared_ptr d_RNG;

public:
    /**\brief  Constructor.
     *\param[in] low  The smallest value the variable can take
     *\param[in] high  The supremum or one more than the max value the variable can take
     *\param[in] r  The generator to use for making the variable
     */
    explicit RandomVariable( type low, type high, RNG::shared_ptr r );

    /**\brief  The casting operator to allow the RandomVariable to be an appropriate
     * rvalue for type.  Everytime it is cast, it generates a new number.
     */
    operator type();
};

/** \brief  Explicit instantiation of a RandomVariable<double>
 * \see RandomVariable
 */
template<>
class RandomVariable<double>
{
public:
    typedef double type;

private:
    type d_Low, d_High;
    RNG::shared_ptr d_RNG;

public:
    explicit RandomVariable( type low, type high, RNG::shared_ptr r );

    operator type();
};

/** \brief  Explicit instantiation of a RandomVariable<float>
 * \see RandomVariable
 */
template<>
class RandomVariable<float>
{
public:
    typedef float type;

private:
    type d_Low, d_High;
    RNG::shared_ptr d_RNG;

public:
    explicit RandomVariable( type low, type high, RNG::shared_ptr r );

    operator type();
};
} // namespace AMP

#include "RNG.inline.h"
#include "RNG.tmpl.h"
#endif
