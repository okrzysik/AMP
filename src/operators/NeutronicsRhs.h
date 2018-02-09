
#ifndef included_AMP_NeutronicsRhs
#define included_AMP_NeutronicsRhs

/* AMP files */
#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "NeutronicsRhsParameters.h"

#include "AMP/utils/InputDatabase.h"

/*Boost files */
#include "AMP/utils/shared_ptr.h"

#include <vector>

namespace AMP {
namespace Operator {

//===========================================================================//
/*!
 * \class NeutronicsRhs
 * \brief A class for representing the neutronics source operator.
 */
//===========================================================================//

class NeutronicsRhs : public Operator
{

public:
    // typedef AMP::LinearAlgebra::VectorVariable<AMP::Mesh::IntegrationPointVariable, 8>
    // HexGaussPointVariable;
    // typedef AMP::shared_ptr<HexGaussPointVariable>      SP_HexGaussPointVariable;
    typedef AMP::shared_ptr<NeutronicsRhsParameters> SP_Parameters;
    typedef AMP::shared_ptr<OperatorParameters> SP_OperatorParameters;
    typedef std::vector<double> Vec_Dbl;
    typedef AMP::shared_ptr<Vec_Dbl> SP_Vec_Dbl;
    typedef AMP::shared_ptr<AMP::Database> SP_Database;

    //! Neutronics Input Types
    enum SourceType { Power, Oxygen, Metal, FissionGas, Isotopes, NUM_SOURCE_TYPES };

private:
public:
    explicit NeutronicsRhs( SP_Parameters parameters );

    /**
     * Empty destructor for NeutronicsRhs
     */
    virtual ~NeutronicsRhs();

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
     * @param f: The result of apply is f = A(u)
     */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /**
      A function to reinitialize this object.
      */
    void reset( const SP_OperatorParameters &parameters ) override;

    // static SP_HexGaussPointVariable createOutputVariable (const std::string & name, int varId =
    // -1);

    void setOutputVariableName( const std::string &name, int varId = -1 );

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override;

    void setTimeStep( int tStep ) { d_timeStep = tStep; }
    void setTimeInSeconds( double seconds );
    void setTimeInDays( double days );

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
    double d_timeStepInSeconds;
    AMP::LinearAlgebra::Variable::shared_ptr d_outputVariable;
    std::vector<Vec_Dbl> d_values;
    double d_secondsPerDay;
    SourceType str2id( const std::string &str );
};
} // namespace Operator
} // namespace AMP

#include "NeutronicsRhs.i.h"

#endif
