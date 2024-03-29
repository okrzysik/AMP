#ifndef included_SolverTestParameters_h
#define included_SolverTestParameters_h

#include "AMP/utils/Database.h"
#include "AMP/utils/UtilityMacros.h"

#include <memory>
#include <string>

namespace AMP::Solver::Test {

struct SolverParameters {

    static std::unique_ptr<AMP::Database> getParameters( const std::string &solver,
                                                         bool use_nested )
    {
        if ( solver == "CG" ) {
            return getCGParameters( use_nested );
        } else if ( solver == "GMRES" ) {
            return getGMRESParameters( use_nested );
        } else if ( solver == "FGMRES" ) {
            return getFGMRESParameters( use_nested );
        } else if ( solver == "BiCGSTAB" ) {
            return getBiCGSTABParameters( use_nested );
        } else if ( solver == "TFQMR" ) {
            return getTFQMRParameters( use_nested );
        } else if ( solver == "PetscFGMRES" ) {
            return getPetscFGMRESParameters( use_nested );
        } else if ( solver == "BoomerAMG" ) {
            return getBoomerAMGParameters( use_nested );
        } else if ( solver == "HyprePCG" ) {
            return getHyprePCGParameters( use_nested );
        } else if ( solver == "ML" ) {
            return getMLParameters( use_nested );
        } else if ( solver == "MueLu" ) {
            return getMueLuParameters( use_nested );
        } else {
            AMP_ERROR( "Unknown solver" );
        }
    }

    static std::unique_ptr<AMP::Database> getCGParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "CGSolver" );
        db->putScalar<bool>( "uses_preconditioner", use_nested );
        db->putScalar<double>( "absolute_tolerance", 1.0e-12 );
        db->putScalar<double>( "relative_tolerance", 1.0e-12 );
        db->putScalar<int>( "print_info_level", 1 );
        db->putScalar<int>( "max_iterations", 100 );
        if ( use_nested )
            db->putScalar<std::string>( "pc_solver_name", "Preconditioner" );
        return db;
    }

    static std::unique_ptr<AMP::Database> getGMRESParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "GMRESSolver" );
        db->putScalar<bool>( "uses_preconditioner", use_nested );
        db->putScalar<double>( "absolute_tolerance", 1.0e-12 );
        db->putScalar<double>( "relative_tolerance", 1.0e-12 );
        db->putScalar<int>( "print_info_level", 1 );
        db->putScalar<int>( "max_iterations", 100 );
        if ( use_nested )
            db->putScalar<std::string>( "pc_solver_name", "Preconditioner" );
        return db;
    }

    static std::unique_ptr<AMP::Database> getFGMRESParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "GMRESSolver" );
        db->putScalar<bool>( "flexible_gmres", true );
        db->putScalar<bool>( "uses_preconditioner", use_nested );
        db->putScalar<double>( "absolute_tolerance", 1.0e-12 );
        db->putScalar<double>( "relative_tolerance", 1.0e-12 );
        db->putScalar<int>( "print_info_level", 1 );
        db->putScalar<int>( "max_iterations", 100 );
        if ( use_nested )
            db->putScalar<std::string>( "pc_solver_name", "Preconditioner" );
        return db;
    }

    static std::unique_ptr<AMP::Database> getBiCGSTABParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "BiCGSTABSolver" );
        db->putScalar<bool>( "uses_preconditioner", use_nested );
        db->putScalar<double>( "absolute_tolerance", 1.0e-12 );
        db->putScalar<double>( "relative_tolerance", 1.0e-12 );
        db->putScalar<int>( "print_info_level", 1 );
        db->putScalar<int>( "max_iterations", 100 );
        if ( use_nested )
            db->putScalar<std::string>( "pc_solver_name", "Preconditioner" );
        return db;
    }

    static std::unique_ptr<AMP::Database> getTFQMRParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "TFQMRSolver" );
        db->putScalar<bool>( "uses_preconditioner", use_nested );
        db->putScalar<double>( "absolute_tolerance", 1.0e-12 );
        db->putScalar<double>( "relative_tolerance", 1.0e-12 );
        db->putScalar<int>( "print_info_level", 1 );
        db->putScalar<int>( "max_iterations", 100 );
        if ( use_nested )
            db->putScalar<std::string>( "pc_solver_name", "Preconditioner" );
        return db;
    }

    static std::unique_ptr<AMP::Database> getPetscFGMRESParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "PetscKrylovSolver" );
        db->putScalar<std::string>( "ksp_type", "fgmres" );
        db->putScalar<std::string>( "pc_type", "shell" );
        db->putScalar<std::string>( "pc_side", "RIGHT" );
        db->putScalar<std::string>( "KSPOptions",
                                    "-ksp_monitor -ksp_converged_reason -ksp_max_it 25 -ksp_rtol "
                                    "1.0e-12 -ksp_atol 1.0e-12" );
        db->putScalar<bool>( "uses_preconditioner", use_nested );
        db->putScalar<double>( "absolute_tolerance", 1.0e-12 );
        db->putScalar<double>( "relative_tolerance", 1.0e-12 );
        db->putScalar<int>( "print_info_level", 1 );
        db->putScalar<int>( "max_iterations", 100 );
        if ( use_nested )
            db->putScalar<std::string>( "pc_solver_name", "Preconditioner" );
        return db;
    }

    static std::unique_ptr<AMP::Database> getHyprePCGParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "HyprePCGSolver" );
        db->putScalar<bool>( "uses_preconditioner", use_nested );
        db->putScalar<double>( "absolute_tolerance", 1.0e-12 );
        db->putScalar<double>( "relative_tolerance", 1.0e-12 );
        db->putScalar<int>( "print_info_level", 1 );
        db->putScalar<int>( "max_iterations", 100 );
        if ( use_nested )
            db->putScalar<std::string>( "pc_solver_name", "Preconditioner" );
        return db;
    }

    static std::unique_ptr<AMP::Database> getBoomerAMGParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "BoomerAMGSolver" );
        db->putScalar<int>( "max_iterations", use_nested ? 1 : 25 );
        db->putScalar<int>( "min_coarse_size", 10 );
        db->putScalar<int>( "relax_type", 16 );
        db->putScalar<int>( "coarsen_type", 10 );
        db->putScalar<int>( "inter_type", 17 );
        db->putScalar<int>( "cycle_type", 1 );
        db->putScalar<int>( "relax_order", 0 );
        db->putScalar<double>( "strong_threshold", 0.5 );
        return db;
    }

    static std::unique_ptr<AMP::Database> getMLParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "TrilinosMLSolver" );
        db->putScalar<int>( "max_iterations", use_nested ? 1 : 25 );
        db->putScalar<bool>( "problem_symmetric", true );
        db->putScalar<std::string>( "smoother_preorpost", "both" );
        db->putScalar<int>( "smoother_sweeps", 2 );
        db->putScalar<int>( "coarse_maxsize", 10 );
        return db;
    }

    static std::unique_ptr<AMP::Database> getMueLuParameters( bool use_nested )
    {
        auto db = std::make_unique<AMP::Database>( "LinearSolver" );
        db->putScalar<std::string>( "name", "TrilinosMueLuSolver" );
        db->putScalar<int>( "max_iterations", use_nested ? 1 : 25 );
        return db;
    }
};

} // namespace AMP::Solver::Test

#endif
