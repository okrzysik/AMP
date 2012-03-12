
#include "utils/Utilities.h"
#include "MLoptions.h"

namespace AMP {
  namespace Solver {

    MLoptions :: MLoptions(const boost::shared_ptr<AMP::Database> &db)
    {
      d_problemType = db->getStringWithDefault("problem_type", "SA");
      addDefaults(d_problemType, db);

      d_maxLevels = db->getInteger("max_levels");
      d_pdeEquations = db->getInteger("PDE_equations");
      d_precType = db->getString("prec_type");

      // increasing sets level 0 to be finest 
      // decreasing sets level (maxLevels - 1) to be the finest.
      d_increasingDecreasing = db->getString("increasingordecreasing");
      d_aggregationType = db->getString("aggregationtype");
      d_aggregationDampingFactor = db->getDouble("aggregation_dampingfactor");
      d_aggregationThreshold = db->getDouble("aggregationthreshold");

      //This is only for METIS and ParMETIS.
      d_nodesPerAggregate = db->getInteger("aggregation_nodes_per_aggregate");
      d_nextLevelAggregatesPerProcess = db->getInteger("aggregation_nextlevel_aggregates_per_process");

      d_eigenAnalysisType = db->getString("eigen-analysis_type");
      d_eigenAnalysisIterations = db->getInteger("eigen-analysis_iterations");
      d_enableEnergyMinimization = db->getBool("energy_minimization_enable");

      d_smootherType = db->getString("smoothertype");
      d_smootherSweeps = db->getInteger("smoother_sweeps");
      d_smootherDampingFactor = db->getDouble("smoother_dampingfactor");
      d_prePost = db->getString("smoother_preorpost");

      d_coarseMaxSize = db->getInteger("coarse_maxsize");
      d_coarseType = db->getString("coarse_type");

      // Should we give ML a list of coordinates?
      d_aggregationAuxEnable = db->getBool("aggregation_aux_enable");
      d_aggregationAuxThreshold = db->getDouble("aggregation_aux_threshold");

      // Should we add vectors to default ML null space
      d_nullSpaceType = db->getString("null_space_type");
      d_nullSpaceDimension = db->getInteger("null_space_dimension");
      d_nullSpaceAddDefaultVectors = db->getBool("null_space_add_default_vectors");
    }

    void MLoptions :: addDefaults(const std::string & problemType, const boost::shared_ptr<AMP::Database> &db) 
    {
      if(problemType == "SA") {
        if(!(db->keyExists("max_levels"))) {
          db->putInteger("max_levels", 10);
        }
        if(!(db->keyExists("PDE_equations"))) {
          db->putInteger("PDE_equations", 1);
        }
        if(!(db->keyExists("prec_type"))) {
          db->putString("prec_type", "MGV");
        }

        if(!(db->keyExists("increasingordecreasing"))) {
          db->putString("increasingordecreasing", "increasing");
        }
        if(!(db->keyExists("aggregationtype"))) {
          db->putString("aggregationtype", "Uncoupled-MIS");
        }
        if(!(db->keyExists("aggregation_dampingfactor"))) {
          db->putDouble("aggregation_dampingfactor", (4.0/3.0));
        }
        if(!(db->keyExists("aggregationthreshold"))) {
          db->putDouble("aggregationthreshold", 0.0);
        }
        //This is only for METIS and ParMETIS.
        if(!(db->keyExists("aggregation_nodes_per_aggregate"))) {
          db->putInteger("aggregation_nodes_per_aggregate", 1);
        }
        if(!(db->keyExists("aggregation_nextlevel_aggregates_per_process"))) {
          db->putInteger("aggregation_nextlevel_aggregates_per_process", 128);
        }

        if(!(db->keyExists("eigen-analysis_type"))) {
          db->putString("eigen-analysis_type", "cg");
        }
        if(!(db->keyExists("eigen-analysis_iterations"))) {
          db->putInteger("eigen-analysis_iterations", 10);
        }
        if(!(db->keyExists("energy_minimization_enable"))) {
          db->putBool("energy_minimization_enable", false);
        }

        if(!(db->keyExists("smoothertype"))) {
          db->putString("smoothertype", "symmetric Gauss-Seidel");
        }
        if(!(db->keyExists("smoother_sweeps"))) {
          db->putInteger("smoother_sweeps", 2);
        }
        if(!(db->keyExists("smoother_dampingfactor"))) {
          db->putDouble("smoother_dampingfactor", 1.0);
        }
        if(!(db->keyExists("smoother_preorpost"))) {
          db->putString("smoother_preorpost", "both");
        }

        if(!(db->keyExists("coarse_maxsize"))) {
          db->putInteger("coarse_maxsize", 128);
        }
        if(!(db->keyExists("coarse_type"))) {
          db->putString("coarse_type", "Amesos-KLU");
        }
        if(!(db->keyExists("aggregation_aux_enable"))) {
          db->putBool("aggregation_aux_enable", false);
        }
        if(!(db->keyExists("aggregation_aux_threshold"))) {
          db->putDouble("aggregation_aux_threshold", 0.0);
        }
        if(!(db->keyExists("null_space_type"))) {
          db->putString("null_space_type", "default vectors");
        }
        if(!(db->keyExists("null_space_dimension"))) {
          db->putInteger("null_space_dimension", 3);
        }
        if(!(db->keyExists("null_space_add_default_vectors"))) {
          db->putBool("null_space_add_default_vectors", true);
        }
      } else if(problemType == "NSSA") {
        if(!(db->keyExists("max_levels"))) {
          db->putInteger("max_levels", 10);
        }
        if(!(db->keyExists("PDE_equations"))) {
          db->putInteger("PDE_equations", 1);
        }
        if(!(db->keyExists("prec_type"))) {
          db->putString("prec_type", "MGW");
        }

        if(!(db->keyExists("increasingordecreasing"))) {
          db->putString("increasingordecreasing", "increasing");
        }
        if(!(db->keyExists("aggregationtype"))) {
          db->putString("aggregationtype", "Uncoupled-MIS");
        }
        if(!(db->keyExists("aggregation_dampingfactor"))) {
          db->putDouble("aggregation_dampingfactor", (4.0/3.0));
        }
        if(!(db->keyExists("aggregationthreshold"))) {
          db->putDouble("aggregationthreshold", 0.0);
        }
        //This is only for METIS and ParMETIS.
        if(!(db->keyExists("aggregation_nodes_per_aggregate"))) {
          db->putInteger("aggregation_nodes_per_aggregate", 1);
        }
        if(!(db->keyExists("aggregation_nextlevel_aggregates_per_process"))) {
          db->putInteger("aggregation_nextlevel_aggregates_per_process", 128);
        }

        if(!(db->keyExists("eigen-analysis_type"))) {
          db->putString("eigen-analysis_type", "power-method");
        }
        if(!(db->keyExists("eigen-analysis_iterations"))) {
          db->putInteger("eigen-analysis_iterations", 20);
        }
        if(!(db->keyExists("energy_minimization_enable"))) {
          db->putBool("energy_minimization_enable", true);
        }

        if(!(db->keyExists("smoothertype"))) {
          db->putString("smoothertype", "symmetric Gauss-Seidel");
        }
        if(!(db->keyExists("smoother_sweeps"))) {
          db->putInteger("smoother_sweeps", 4);
        }
        if(!(db->keyExists("smoother_dampingfactor"))) {
          db->putDouble("smoother_dampingfactor", 0.67);
        }
        if(!(db->keyExists("smoother_preorpost"))) {
          db->putString("smoother_preorpost", "post");
        }

        if(!(db->keyExists("coarse_maxsize"))) {
          db->putInteger("coarse_maxsize", 256);
        }
        if(!(db->keyExists("coarse_type"))) {
          db->putString("coarse_type", "Amesos-KLU");
        }
        if(!(db->keyExists("aggregation_aux_enable"))) {
          db->putBool("aggregation_aux_enable", false);
        }
        if(!(db->keyExists("aggregation_aux_threshold"))) {
          db->putDouble("aggregation_aux_threshold", 0.0);
        }
        if(!(db->keyExists("null_space_type"))) {
          db->putString("null_space_type", "default vectors");
        }
        if(!(db->keyExists("null_space_dimension"))) {
          db->putInteger("null_space_dimension", 3);
        }
        if(!(db->keyExists("null_space_add_default_vectors"))) {
          db->putBool("null_space_add_default_vectors", true);
        }
      } else {
        AMP_ERROR("The option, problem_type = \"" << problemType << "\" , is not supported.");
      }
    }

  }
}

