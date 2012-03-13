
#include "utils/Utilities.h"
#include "TrilinosMLSolver.h"
#include "vectors/trilinos/EpetraVector.h"
#include "matrices/Matrix.h"
#include "matrices/trilinos/EpetraMatrix.h"
#include "operators/LinearOperator.h"
#include "operators/TrilinosMatrixShellOperator.h"

namespace AMP {
  namespace Solver {

    TrilinosMLSolver::TrilinosMLSolver()
    {
      d_ml = NULL;
      d_mlAggregate = NULL;
      d_bCreationPhase = true; 
    }

    TrilinosMLSolver::TrilinosMLSolver(boost::shared_ptr<SolverStrategyParameters> parameters):SolverStrategy(parameters)
    {
      d_ml = NULL;
      d_mlAggregate = NULL;
      assert(parameters.get()!=NULL);
      initialize(parameters);
    }

    TrilinosMLSolver::~TrilinosMLSolver()
    {
      if(d_mlAggregate) {
        ML_Aggregate_Destroy(&d_mlAggregate);
        d_mlAggregate = NULL;
      }
      if(d_ml) {
        ML_Destroy(&d_ml);
        d_ml = NULL;
      }
    }

    void
      TrilinosMLSolver::getFromInput(const boost::shared_ptr<AMP::Database> &db)
      {
        d_bUseEpetra = db->getBoolWithDefault("USE_EPETRA", true);

        d_mlOptions.reset(new MLoptions(db));

        if(d_bUseEpetra) {
          convertMLoptionsToTeuchosParameterList();
        }
      }

    void 
      TrilinosMLSolver::convertMLoptionsToTeuchosParameterList() 
      {
        // output level, 0 being silent and 10 verbose
        d_MLParameterList.set("ML output", d_iDebugPrintInfoLevel);

        // maximum number of levels
        d_MLParameterList.set("max levels", d_mlOptions->d_maxLevels);
        d_MLParameterList.set("prec type", d_mlOptions->d_precType);
        d_MLParameterList.set("PDE equations", d_mlOptions->d_pdeEquations);
        d_MLParameterList.set("cycle applications", d_iMaxIterations);

        d_MLParameterList.set("increasing or decreasing", d_mlOptions->d_increasingDecreasing);
        d_MLParameterList.set("aggregation: type", d_mlOptions->d_aggregationType);
        d_MLParameterList.set("aggregation: damping factor", d_mlOptions->d_aggregationDampingFactor);
        d_MLParameterList.set("aggregation: threshold", d_mlOptions->d_aggregationThreshold);
        d_MLParameterList.set("aggregation: nodes per aggregate", d_mlOptions->d_nodesPerAggregate);
        d_MLParameterList.set("aggregation: next-level aggregates per process", d_mlOptions->d_nextLevelAggregatesPerProcess);

        d_MLParameterList.set("eigen-analysis: type", d_mlOptions->d_eigenAnalysisType);
        d_MLParameterList.set("eigen-analysis: iterations", d_mlOptions->d_eigenAnalysisIterations);

        d_MLParameterList.set("smoother: sweeps", d_mlOptions->d_smootherSweeps);
        d_MLParameterList.set("smoother: damping factor", d_mlOptions->d_smootherDampingFactor);
        d_MLParameterList.set("smoother: pre or post", d_mlOptions->d_prePost);
        d_MLParameterList.set("smoother: type", d_mlOptions->d_smootherType);

        d_MLParameterList.set("energy minimization: enable", d_mlOptions->d_enableEnergyMinimization);

        d_MLParameterList.set("coarse: type", d_mlOptions->d_coarseType);
        d_MLParameterList.set("coarse: max size", d_mlOptions->d_coarseMaxSize);

        d_MLParameterList.set("aggregation: aux: enable", d_mlOptions->d_aggregationAuxEnable);
        d_MLParameterList.set("aggregation: aux: threshold", d_mlOptions->d_aggregationAuxThreshold);

        d_MLParameterList.set("null space: type", d_mlOptions->d_nullSpaceType);
        d_MLParameterList.set("null space: dimension", d_mlOptions->d_nullSpaceDimension);
        d_MLParameterList.set("null space: add default vectors", d_mlOptions->d_nullSpaceAddDefaultVectors);
      }

    void
      TrilinosMLSolver::initialize(boost::shared_ptr<SolverStrategyParameters> const parameters)
      {
        getFromInput(parameters->d_db);

        if(d_pOperator.get() != NULL)
        {
          registerOperator(d_pOperator);
        }

      }

    void
      TrilinosMLSolver::registerOperator(const boost::shared_ptr<AMP::Operator::Operator> op)
      {
        d_pOperator = op;
        AMP_INSIST(d_pOperator.get()!=NULL,"ERROR: TrilinosMLSolver::initialize() operator cannot be NULL");

        if(d_bUseEpetra) {
          // Compute coordinates to give to ML if requested
          if( d_mlOptions->d_aggregationAuxEnable || 
              d_mlOptions->d_nullSpaceType == "from coordinates")
          {
            computeCoordinates( op );
            AMP_ASSERT( d_x_values.size() != 0 );
            AMP_ASSERT( d_y_values.size() != 0 );
            AMP_ASSERT( d_z_values.size() != 0 );
            d_MLParameterList.set("x-coordinates", &d_x_values[0]);
            d_MLParameterList.set("y-coordinates", &d_y_values[0]);
            d_MLParameterList.set("z-coordinates", &d_z_values[0]);
          }

          // Compute null space manually if requested
          if( d_mlOptions->d_nullSpaceType == "pre-computed" )
          {
            computeNullSpace( op );
            AMP_ASSERT( d_null_space.size() != 0 );
            d_MLParameterList.set("null space: vectors",&d_null_space[0]);
          }

          AMP_INSIST( d_mlOptions->d_nullSpaceType == "from coordinates" ?
                      d_mlOptions->d_pdeEquations  == 3 : true,
                      "Null space construction only available for mechanics (PDE_equations=3)");
          AMP_INSIST( d_mlOptions->d_nullSpaceType == "pre-computed" ?
                      d_mlOptions->d_pdeEquations  == 3 : true,
                      "Null space construction only available for mechanics (PDE_equations=3)");

          boost::shared_ptr<AMP::Operator::LinearOperator> linearOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>(d_pOperator);
          AMP_INSIST(linearOperator.get() != NULL, "linearOperator cannot be NULL");

          boost::shared_ptr<AMP::LinearAlgebra::EpetraMatrix> pMatrix = boost::dynamic_pointer_cast<
            AMP::LinearAlgebra::EpetraMatrix>(linearOperator->getMatrix());
          AMP_INSIST(pMatrix.get()!=NULL, "pMatrix cannot be NULL");

          d_mlSolver.reset( new ML_Epetra::MultiLevelPreconditioner(pMatrix->getEpetra_CrsMatrix(), d_MLParameterList, false));


        } else {
          boost::shared_ptr<AMP::Operator::TrilinosMatrixShellOperator> matShellOperator = boost::dynamic_pointer_cast<
            AMP::Operator::TrilinosMatrixShellOperator>(d_pOperator);
          AMP_ASSERT(matShellOperator.get() != NULL);

          size_t matSize = matShellOperator->getMatrixSize();
          ML_Create(&d_ml, d_mlOptions->d_maxLevels);

          if((d_mlOptions->d_increasingDecreasing) == "increasing") {
            ML_Init_Amatrix(d_ml, 0, matSize, matSize, d_pOperator.get());
            ML_Set_Amatrix_Getrow(d_ml, 0, &(AMP::Operator::TrilinosMatrixShellOperator::getRow), NULL, matSize);
            ML_Set_Amatrix_Matvec(d_ml, 0, &(AMP::Operator::TrilinosMatrixShellOperator::matVec));
          } else {
            AMP_ERROR("The option, increasingordecreasing = \"" << (d_mlOptions->d_increasingDecreasing) << "\" , is not supported.");
          }
        }

        d_bCreationPhase = true;
      }

    void
      TrilinosMLSolver::resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params)
      {
        AMP_INSIST((d_pOperator.get() != NULL), "ERROR: TrilinosMLSolver::resetOperator() operator cannot be NULL");

        d_pOperator->reset(params);

        reset( boost::shared_ptr<SolverStrategyParameters>() );
      }

    void
      TrilinosMLSolver::reset(boost::shared_ptr<SolverStrategyParameters> )
      {
        if(!d_bCreationPhase)
        {
          if(d_bUseEpetra) {
            d_mlSolver->DestroyPreconditioner();
          } else {
            ML_Aggregate_Destroy(&d_mlAggregate);
            d_mlAggregate = NULL;
          }
        }
        d_bCreationPhase = true;
      }

    void
      TrilinosMLSolver::solve(boost::shared_ptr<AMP::LinearAlgebra::Vector>  f,
          boost::shared_ptr<AMP::LinearAlgebra::Vector>  u)
      {
        // in this case we make the assumption we can access a EpetraMat for now
        AMP_INSIST(d_pOperator.get()!=NULL,"ERROR: TrilinosMLSolver::solve() operator cannot be NULL");

        if(d_bUseZeroInitialGuess) {
          u->zero();
        }

        if(d_bCreationPhase) {
          if(d_bUseEpetra) {
            d_mlSolver->ComputePreconditioner();
            if(d_iDebugPrintInfoLevel > 2)
            {
              d_mlSolver->PrintUnused();
            }
          } else {
            buildML();
          }
          d_bCreationPhase = false;
        }

        boost::shared_ptr <AMP::LinearAlgebra::Vector> r;  

        if( d_iDebugPrintInfoLevel > 1)
        {
          r = f->cloneVector();
          d_pOperator->apply(f,u,r);
          AMP::pout << "TrilinosMLSolver::solve(), L2 norm of residual before solve " <<std::setprecision(15)<< r->L2Norm() << std::endl;
        }

        if( d_iDebugPrintInfoLevel > 2) 
        {
          double solution_norm = u->L2Norm();
          AMP::pout << "TrilinosMLSolver : before solve solution norm: " <<std::setprecision(15)<< solution_norm << std::endl;
        }

        if(d_bUseEpetra) {
          // These functions throw exceptions if this cannot be performed.
          Epetra_Vector &fVec = (AMP::LinearAlgebra::EpetraVector::view ( f ))->castTo<AMP::LinearAlgebra::EpetraVector>().getEpetra_Vector();
          Epetra_Vector &uVec = (AMP::LinearAlgebra::EpetraVector::view ( u ))->castTo<AMP::LinearAlgebra::EpetraVector>().getEpetra_Vector();

          d_mlSolver->ApplyInverse(fVec, uVec);
        } else {
          double * uArr = u->getRawDataBlock<double>();
          double * fArr = f->getRawDataBlock<double>();

          ML_Iterate(d_ml, uArr, fArr);
        }

        // we are forced to update the state of u here
        // as Epetra is not going to change the state of a managed vector
        // an example where this will and has caused problems is when the
        // vector is a petsc managed vector being passed back to PETSc
        if ( u->isA<AMP::LinearAlgebra::DataChangeFirer>() )
        {
          u->castTo<AMP::LinearAlgebra::DataChangeFirer>().fireDataChange();
        }

        if( d_iDebugPrintInfoLevel > 2) 
        {
          double solution_norm = u->L2Norm();
          AMP::pout << "TrilinosMLSolver : after solve solution norm: " <<std::setprecision(15)<< solution_norm << std::endl;
        }

        // for debugging in extreme cases, clone a temporary vector and compute the residual
        if( d_iDebugPrintInfoLevel > 1)
        {
          d_pOperator->apply(f,u,r);
          AMP::pout << "TrilinosMLSolver::solve(), L2 norm of residual after solve " <<std::setprecision(15)<< r->L2Norm() << std::endl;    
        }

      }

    void 
      TrilinosMLSolver::buildML()
      {
        ML_Set_MaxIterations(d_ml, d_iMaxIterations);
        ML_Set_PrintLevel(d_iDebugPrintInfoLevel);
        ML_Set_OutputLevel(d_ml, d_iDebugPrintInfoLevel);
        if(d_iDebugPrintInfoLevel) {
          ML_Set_ResidualOutputFrequency(d_ml, 1);
        }

        ML_Aggregate_Create(&d_mlAggregate);

        d_mlAggregate->num_PDE_eqns = d_mlOptions->d_pdeEquations;
        d_mlAggregate->nullspace_dim = d_mlOptions->d_pdeEquations;

        ML_Aggregate_Set_MaxCoarseSize(d_mlAggregate, (d_mlOptions->d_coarseMaxSize));
        if((d_mlOptions->d_aggregationType) == "Uncoupled-MIS") {
          ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(d_mlAggregate);
        } else {
          AMP_ERROR("The option, aggregationtype = \"" << (d_mlOptions->d_aggregationType) << "\" , is not supported.");
        }

        int nlevels = ML_Gen_MGHierarchy_UsingAggregation(d_ml, 0, ML_INCREASING, d_mlAggregate);
        AMP::pout<<"Number of actual levels : "<< nlevels <<std::endl;

        if((d_mlOptions->d_smootherType) == "symmetric Gauss-Seidel") {
          for(int lev = 0; lev < (nlevels - 1); lev++) {
            if((d_mlOptions->d_prePost) == "pre") {
              ML_Gen_Smoother_SymGaussSeidel(d_ml, lev, ML_PRESMOOTHER,
                  (d_mlOptions->d_smootherSweeps), (d_mlOptions->d_smootherDampingFactor));
            } else if((d_mlOptions->d_prePost) == "post") {
              ML_Gen_Smoother_SymGaussSeidel(d_ml, lev, ML_POSTSMOOTHER,
                  (d_mlOptions->d_smootherSweeps), (d_mlOptions->d_smootherDampingFactor));
            } else if((d_mlOptions->d_prePost) == "both") {
              ML_Gen_Smoother_SymGaussSeidel(d_ml, lev, ML_BOTH,
                  (d_mlOptions->d_smootherSweeps), (d_mlOptions->d_smootherDampingFactor));
            } else {
              AMP_ERROR("The option, smoother_preorpost = \"" << (d_mlOptions->d_prePost) << "\" , is not supported.");
            }
          }
        } else {
          AMP_ERROR("The option, smoothertype = \"" << (d_mlOptions->d_smootherType) << "\" , is not supported.");
        }

        if((d_mlOptions->d_coarseType) == "Amesos-KLU") {
          ML_Gen_Smoother_Amesos(d_ml, (nlevels - 1), ML_AMESOS_KLU, -1, 0.0);
        } else {
          AMP_ERROR("The option, coarse_type = \"" << (d_mlOptions->d_coarseType) << "\" , is not supported.");
        }

        if((d_mlOptions->d_precType) == "MGV") {
          ML_Gen_Solver(d_ml, ML_MGV, 0, (nlevels - 1));
        } else if((d_mlOptions->d_precType) == "MGW") {
          ML_Gen_Solver(d_ml, ML_MGW, 0, (nlevels - 1));
        } else {
          AMP_ERROR("The option, prec_type = \"" << (d_mlOptions->d_precType) << "\" , is not supported.");
        }
      }

    void 
        TrilinosMLSolver::computeCoordinates( const boost::shared_ptr<AMP::Operator::Operator> op )
        {
            // Get mesh adapter for this operator
            AMP::Mesh::MeshManager::Adapter::shared_ptr myMesh = op->getMeshAdapter();

            // Resize vectors to hold node values
            int numNodes = myMesh->numLocalNodes();
            d_x_values.resize(numNodes,0.0);
            d_y_values.resize(numNodes,0.0);
            d_z_values.resize(numNodes,0.0);

            // Get node iterators
            AMP::Mesh::MeshManager::Adapter::NodeIterator thisNode = myMesh->beginNode();
            AMP::Mesh::MeshManager::Adapter::NodeIterator  endNode = myMesh->endNode();

            int nodeCounter = 0;
            for( ; thisNode != endNode; ++thisNode )
            {
                d_x_values[nodeCounter] = (*thisNode).x();
                d_y_values[nodeCounter] = (*thisNode).y();
                d_z_values[nodeCounter] = (*thisNode).z();
                nodeCounter++;
            }
        }

    void 
        TrilinosMLSolver::computeNullSpace( const boost::shared_ptr<AMP::Operator::Operator> op )
        {
            // Get mesh adapter for this operator
            AMP::Mesh::MeshManager::Adapter::shared_ptr myMesh = op->getMeshAdapter();
            int numPDE = d_mlOptions->d_pdeEquations;
            int dimNS  = d_mlOptions->d_nullSpaceDimension;


            // Resize vectors to hold node values
            int numNodes = myMesh->numLocalNodes();
            int vecLength = numPDE*numNodes;
            d_null_space.resize(dimNS*vecLength,0.0);

            // Get node iterators
            AMP::Mesh::MeshManager::Adapter::NodeIterator thisNode = myMesh->beginNode();
            AMP::Mesh::MeshManager::Adapter::NodeIterator  endNode = myMesh->endNode();

            double thisX;
            double thisY;
            double thisZ;
            int nodeCounter = 0;
            int offset = 0;
            int vecOffset = 0;
            for( ; thisNode != endNode; ++thisNode )
            {
                thisX = (*thisNode).x();
                thisY = (*thisNode).y();
                thisZ = (*thisNode).z();

                int dof = numPDE * nodeCounter;

                // Constant vector for each PDE
                for( int i=0; i<numPDE; ++i )
                {
                    vecOffset = i * vecLength;
                    offset    = vecOffset + dof + i;
                    d_null_space[offset] = 1.0;
                }

                // Rotation around X
                vecOffset = 3*vecLength;
                offset    = vecOffset + dof;
                d_null_space[offset + 1] = -thisZ;
                d_null_space[offset + 2] =  thisY;

                // Rotation around Y
                vecOffset = 4*vecLength;
                offset    = vecOffset + dof;
                d_null_space[offset]     =  thisZ;
                d_null_space[offset + 2] = -thisX;

                // Rotation around Z
                vecOffset = 5*vecLength;
                offset    = vecOffset + dof;
                d_null_space[offset]     = -thisY;
                d_null_space[offset + 1] =  thisX;

                nodeCounter++;
            }
        }
  }
}


