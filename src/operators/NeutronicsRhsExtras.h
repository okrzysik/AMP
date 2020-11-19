#ifndef included_AMP_NeutronicsRhsExtras
#define included_AMP_NeutronicsRhsExtras

#include "AMP/operators/NeutronicsRhsExtrasParameters.h"
#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include <memory>

#include <vector>

namespace AMP {
namespace Operator {

//===========================================================================//
/*!
 * \class NeutronicsRhsExtras
 * \brief A class for representing the neutronics source operator.
 */
//===========================================================================//

class NeutronicsRhsExtras : public Operator
{

public:
    typedef std::shared_ptr<NeutronicsRhsExtrasParameters> SP_Parameters;
    typedef std::shared_ptr<OperatorParameters> SP_OperatorParameters;
    typedef std::vector<double> Vec_Dbl;
    typedef std::shared_ptr<Vec_Dbl> SP_Vec_Dbl;
    typedef std::shared_ptr<AMP::Database> SP_Database;
    typedef std::vector<double> Vec_Dbl1;
    typedef std::vector<Vec_Dbl1> Vec_Dbl2;
    typedef std::vector<Vec_Dbl2> Vec_Dbl3;

    //! Neutronics Input Types
    enum SourceType { Isotopes, Elements, NUM_SOURCE_TYPES };

private:
public:
    explicit NeutronicsRhsExtras( SP_Parameters parameters );

    /**
     * Empty destructor for NeutronicsRhsExtras
     */
    virtual ~NeutronicsRhsExtras();

    //! Return the name of the operator
    std::string type() const override { return "NeutronicsRhsExtras"; }

    /**
     * Print out all members of integrator instance to given output stream.
     */
    void printClassData( std::ostream &os ) const;

    /**
     * Write out state of object to given database.
     *
     * When assertion checking is active, the database pointer must be non-null.
     */
    void putToDatabase( SP_Database db );

    /**
      The function that computes the residual.
     * @param u: multivector of the state.
     * @param f: The result of apply is r = A(u)
     */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /**
      A function to reinitialize this object.
      */
    void reset( const SP_OperatorParameters &parameters ) override;

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override
    {
        return d_outputVariable;
    }

    void setTimeStep( int tStep ) { d_timeStep = tStep; }
    void setExtrasId( int extrasId ) { d_extrasId = extrasId; }
    void setTimeInSeconds( double seconds );

protected:
    /*
     * Read input data from specified database and initialize class members.
     * If run is from restart, a subset of the restart values may be replaced
     * with those read from input.
     *
     * When assertion checking is active, the database pointer must be non-null.
     */
    void getFromInput( SP_Database db );

    SP_Database d_db;
    bool d_useFixedValue;
    int d_numTimeSteps;
    Vec_Dbl d_timeStepsInDays;
    SourceType d_type;
    Vec_Dbl d_fixedValues;
    int d_timeStep;
    int d_extrasId;
    int d_numExtras;
    std::vector<std::string> d_extrasName;
    double d_timeStepInSeconds;
    AMP::LinearAlgebra::Variable::shared_ptr d_outputVariable;
    Vec_Dbl3 d_values;
    double d_secondsPerDay;
    SourceType str2id( const std::string &str );
};
} // namespace Operator
} // namespace AMP

#include "NeutronicsRhsExtras.i.h"

#endif
