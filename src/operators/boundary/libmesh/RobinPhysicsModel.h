#ifndef included_AMP_RobinPhysicsModel
#define included_AMP_RobinPhysicsModel

#include <iostream>
#include <string>
#include <vector>

#include "materials/Material.h"
#include "materials/Property.h"
#include "operators/ElementPhysicsModel.h"
#include "utils/shared_ptr.h"
//#include "RobinPhysicsModelParameters.h"

namespace AMP {
namespace Operator {

typedef ElementPhysicsModelParameters RobinPhysicsModelParameters;

class RobinPhysicsModel : public ElementPhysicsModel
{
public:
    explicit RobinPhysicsModel( const AMP::shared_ptr<RobinPhysicsModelParameters> &params )
        : ElementPhysicsModel( params ), d_numActiveVariables( 0 )
    {
        reset( params );
    }

    virtual ~RobinPhysicsModel() {}

    virtual void reset( const AMP::shared_ptr<RobinPhysicsModelParameters> &params )
    {
        if ( params->d_db->keyExists( "Number_Active_Variables" ) ) {
            d_numActiveVariables = ( params->d_db )->getInteger( "Number_Active_Variables" );
        }
        AMP::shared_ptr<AMP::Database> activeDb;
        if ( params->d_db->keyExists( "ActiveInputVariables" ) ) {
            activeDb               = params->d_db->getDatabase( "ActiveInputVariables" );
            unsigned int numactive = activeDb->getAllKeys().size();
            if ( params->d_db->keyExists( "Number_Active_Variables" ) ) {
                AMP_INSIST( numactive == d_numActiveVariables,
                            "Number of active variables disagrees with Number_Active_Variables." );
            }
            d_numActiveVariables = numactive;
        }

        d_activeVariableNames.resize( d_numActiveVariables );
        for ( unsigned int var = 0; var < d_numActiveVariables; var++ ) {
            char key[100];
            sprintf( key, "ActiveVariable_%u", var );
            std::string varName        = activeDb->getString( key );
            d_activeVariableNames[var] = varName;
        }
    }

    virtual void getConductance( std::vector<double> &beta,
                                 std::vector<double> &gamma,
                                 const std::vector<std::vector<double>> &arguments ) = 0;

    const std::vector<std::string> getVariableName() { return d_activeVariableNames; }

protected:
    unsigned int d_numActiveVariables; /**< Number of Active Variables. */

    std::vector<std::string>
        d_activeVariableNames; /**< A list of strings to store Active Variable names. */

private:
};
}
}

#endif
