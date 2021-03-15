
#include "MLoptions.h"
#include "AMP/utils/Utilities.h"

namespace AMP {
namespace Solver {

MLoptions::MLoptions( std::shared_ptr<AMP::Database> db )
    : d_problemType( db->getWithDefault<std::string>( "problem_type", "SA" ) )
{
    addDefaults( d_problemType, db );

    d_maxLevels    = db->getScalar<int>( "max_levels" );
    d_pdeEquations = db->getScalar<int>( "PDE_equations" );
    d_precType     = db->getString( "prec_type" );

    // increasing sets level 0 to be finest
    // decreasing sets level (maxLevels - 1) to be the finest.
    d_increasingDecreasing     = db->getString( "increasingordecreasing" );
    d_aggregationType          = db->getString( "aggregationtype" );
    d_aggregationDampingFactor = db->getScalar<double>( "aggregation_dampingfactor" );
    d_aggregationThreshold     = db->getScalar<double>( "aggregationthreshold" );

    // This is only for METIS and ParMETIS.
    d_nodesPerAggregate = db->getScalar<int>( "aggregation_nodes_per_aggregate" );
    d_nextLevelAggregatesPerProcess =
        db->getScalar<int>( "aggregation_nextlevel_aggregates_per_process" );

    d_eigenAnalysisType        = db->getString( "eigen-analysis_type" );
    d_eigenAnalysisIterations  = db->getScalar<int>( "eigen-analysis_iterations" );
    d_enableEnergyMinimization = db->getScalar<bool>( "energy_minimization_enable" );

    d_smootherType          = db->getString( "smoothertype" );
    d_smootherSweeps        = db->getScalar<int>( "smoother_sweeps" );
    d_smootherDampingFactor = db->getScalar<double>( "smoother_dampingfactor" );
    d_prePost               = db->getString( "smoother_preorpost" );

    d_coarseMaxSize = db->getScalar<int>( "coarse_maxsize" );
    d_coarseType    = db->getString( "coarse_type" );

    // Should we give ML a list of coordinates?
    d_aggregationAuxEnable    = db->getScalar<bool>( "aggregation_aux_enable" );
    d_aggregationAuxThreshold = db->getScalar<double>( "aggregation_aux_threshold" );

    // Should we add vectors to default ML null space
    d_nullSpaceType              = db->getString( "null_space_type" );
    d_nullSpaceDimension         = db->getScalar<int>( "null_space_dimension" );
    d_nullSpaceAddDefaultVectors = db->getScalar<bool>( "null_space_add_default_vectors" );
}

template<class TYPE>
static inline void addEntry( std::shared_ptr<AMP::Database> db, std::string_view key, TYPE value )
{
    if ( !db->keyExists( key ) )
        db->putScalar( key, value );
}
void MLoptions::addDefaults( const std::string &problemType, std::shared_ptr<AMP::Database> db )
{
    if ( problemType == "SA" ) {
        addEntry( db, "max_levels", 10 );
        addEntry( db, "PDE_equations", 1 );
        addEntry( db, "prec_type", "MGV" );
        addEntry( db, "increasingordecreasing", "increasing" );
        addEntry( db, "aggregationtype", "Uncoupled-MIS" );
        addEntry( db, "aggregation_dampingfactor", 4.0 / 3.0 );
        addEntry( db, "aggregationthreshold", 0.0 );
        // This is only for METIS and ParMETIS.
        addEntry( db, "aggregation_nodes_per_aggregate", 1 );
        addEntry( db, "aggregation_nextlevel_aggregates_per_process", 128 );
        addEntry( db, "eigen-analysis_type", "cg" );
        addEntry( db, "eigen-analysis_iterations", 10 );
        addEntry( db, "energy_minimization_enable", false );
        addEntry( db, "smoothertype", "symmetric Gauss-Seidel" );
        addEntry( db, "smoother_sweeps", 2 );
        addEntry( db, "smoother_dampingfactor", 1.0 );
        addEntry( db, "smoother_preorpost", "both" );
        addEntry( db, "coarse_maxsize", 128 );
        addEntry( db, "coarse_type", "Amesos-KLU" );
        addEntry( db, "aggregation_aux_enable", false );
        addEntry( db, "aggregation_aux_threshold", 0.0 );
        addEntry( db, "null_space_type", "default vectors" );
        addEntry( db, "null_space_dimension", 3 );
        addEntry( db, "null_space_add_default_vectors", false );
    } else if ( problemType == "NSSA" ) {
        addEntry( db, "max_levels", 10 );
        addEntry( db, "PDE_equations", 1 );
        addEntry( db, "prec_type", "MGW" );
        addEntry( db, "increasingordecreasing", "increasing" );
        addEntry( db, "aggregationtype", "Uncoupled-MIS" );
        addEntry( db, "aggregation_dampingfactor", ( 4.0 / 3.0 ) );
        addEntry( db, "aggregationthreshold", 0.0 );
        // This is only for METIS and ParMETIS.
        addEntry( db, "aggregation_nodes_per_aggregate", 1 );
        addEntry( db, "aggregation_nextlevel_aggregates_per_process", 128 );
        addEntry( db, "eigen-analysis_type", "power-method" );
        addEntry( db, "eigen-analysis_iterations", 20 );
        addEntry( db, "energy_minimization_enable", true );
        addEntry( db, "smoothertype", "symmetric Gauss-Seidel" );
        addEntry( db, "smoother_sweeps", 4 );
        addEntry( db, "smoother_dampingfactor", 0.67 );
        addEntry( db, "smoother_preorpost", "post" );
        addEntry( db, "coarse_maxsize", 256 );
        addEntry( db, "coarse_type", "Amesos-KLU" );
        addEntry( db, "aggregation_aux_enable", false );
        addEntry( db, "aggregation_aux_threshold", 0.0 );
        addEntry( db, "null_space_type", "default vectors" );
        addEntry( db, "null_space_dimension", 3 );
        addEntry( db, "null_space_add_default_vectors", false );
    } else {
        AMP_ERROR( "The option, problem_type = \"" << problemType << "\" , is not supported." );
    }
}
} // namespace Solver
} // namespace AMP
