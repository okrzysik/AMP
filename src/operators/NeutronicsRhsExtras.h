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

namespace AMP::Operator {

//===========================================================================//
/*!
 * \class NeutronicsRhsExtras
 * \brief A class for representing the neutronics source operator.
 */
//===========================================================================//

class NeutronicsRhsExtras : public Operator
{
public:
    //! Neutronics Input Types
    enum SourceType { Isotopes, Elements, NUM_SOURCE_TYPES };

public:
    explicit NeutronicsRhsExtras( std::shared_ptr<NeutronicsRhsExtrasParameters> parameters );

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
    void putToDatabase( std::shared_ptr<AMP::Database> db );

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
    void reset( std::shared_ptr<const OperatorParameters> parameters ) override;

    std::shared_ptr<AMP::LinearAlgebra::Variable> getOutputVariable() override
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
    void getFromInput( std::shared_ptr<AMP::Database> db );

    std::shared_ptr<AMP::Database> d_db;
    bool d_useFixedValue;
    int d_numTimeSteps;
    std::vector<double> d_timeStepsInDays;
    SourceType d_type;
    std::vector<double> d_fixedValues;
    int d_timeStep;
    int d_extrasId;
    int d_numExtras;
    std::vector<std::string> d_extrasName;
    double d_timeStepInSeconds;
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_outputVariable;
    std::vector<std::vector<std::vector<double>>> d_values;
    double d_secondsPerDay;
    SourceType str2id( const std::string &str );
};
} // namespace AMP::Operator

#endif
