#include "MoertelOperatorBuilder.h"

#include "EpetraMatrixOperator.h"
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"
#include "vectors/trilinos/EpetraVectorEngine.h"
#include "vectors/trilinos/ManagedEpetraVector.h"


namespace AMP {
namespace Operator {
/*
  void MoertelOperatorBuilder::transposeMatrix ( const Epetra_CrsMatrix *in , Epetra_CrsMatrix *&out )
  {
    std::map<int , std::map<int , double> >  aij;
    std::vector<int>  sizes ( in->NumGlobalCols());
    for ( int i = 0 ; i != in->NumMyRows() ; i++ )
    {
      int numItems;
      double *vals;
      int * cols;
      in->ExtractMyRowView ( i , numItems , vals , cols );
      for ( int j = 0 ; j != numItems ; j++ )
      {
        int realCol = in->ColMap().GID ( cols[j] );
        aij[realCol][i] = vals[j];
        sizes[realCol]++;
      }
    }

    out = new Epetra_CrsMatrix ( Copy , in->DomainMap() , &(sizes[0]) );
    std::map<int , std::map<int , double> >::iterator cur_row = aij.begin();
    while ( cur_row != aij.end() )
    {
      std::vector<double>  valsIn; 
      std::vector<int>     colsIn; 
      valsIn.reserve ( cur_row->second.size() );
      colsIn.reserve ( cur_row->second.size() );
      std::map<int,double>::iterator cur_aij = cur_row->second.begin();
      while ( cur_aij != cur_row->second.end() )
      {
        colsIn.push_back ( cur_aij->first );
        valsIn.push_back ( cur_aij->second );
        cur_aij++;
      }
      int t1 = cur_row->first;
      int t2 = colsIn.size();
      out->InsertGlobalValues ( t1 , t2 , &(valsIn[0]) , &(colsIn[0]) );
      cur_row++;
    }
    out->FillComplete ( in->RangeMap() , in->DomainMap() );
  }

  void MoertelOperatorBuilder::subMatrix ( const Epetra_CrsMatrix *in , Epetra_CrsMatrix *&out , int row_start , int row_end , int col_start , int col_end )
  {
    int  numRows = row_end - row_start;
    std::vector<double *>  valsView ( numRows );
    std::vector<int *>     colsView ( numRows );
    std::vector<std::vector<int> >      cols ( numRows );
    std::vector<int>       numCols ( numRows );
    std::vector<int>       numColsView ( numRows );
    std::vector<int>       colOffs ( numRows );
    int maxid = 0;
    for ( int i = row_start ; i != row_end ; i++ )
    {
      int localID = i - row_start;
      in->ExtractMyRowView ( i , numColsView[localID] , valsView[localID] , colsView[localID] );
      colOffs[localID] = 0;
      numCols[localID] = 0;
      int numCurCols = numColsView[localID];
      cols[localID].resize ( numCurCols );
      for ( int j = 0 ; j != numCurCols ; j++ )
      {
        int curColID = colsView[localID][j];
        int oldId = in->ColMap().GID ( curColID );
        int newId = oldId - col_start;
        maxid = std::max ( newId , maxid );
        cols[localID][j] = newId;
        if ( oldId < col_start )
        {
          colOffs[localID]++;
        }
        else if ( oldId < col_end )
        {
          numCols[localID]++;
        }
      }
    }
    Epetra_Map newRangeMap ( row_end - row_start , 0 , d_Comm );
    Epetra_Map newDomainMap ( col_end - col_start , 0 , d_Comm );
    out = new Epetra_CrsMatrix ( Copy , newRangeMap , &(numCols[0]) );
    for ( int i = 0 ; i != numRows ; i++ )
    {
      out->InsertGlobalValues ( i , numCols[i] , valsView[i] + colOffs[i] , &(cols[i][0]) + colOffs[i] );
    }
    out->FillComplete ( newDomainMap , newRangeMap );
  }

  void MoertelOperatorBuilder::getNodes ( AMP::Mesh::MeshManager::Adapter::shared_ptr mesh , NodeIterator begin , NodeIterator end , MOERTEL::Interface &interface , AMP::Mesh::DOFMap::shared_ptr  d_map , int side , size_t offset )
  {
    while ( begin != end )
    {
      double pos[3];
      pos[0] = begin->x(); pos[1] = begin->y(); pos[2] = begin->z();
      int dof[1];
      dof[0] = d_map->getGlobalID ( begin->globalID() , 0 ) + offset;
      bool onBoundary = mesh->getBoundaryIds ( *begin ).size() > 1;
      d_NodeRenumber [ begin->globalID()+offset ] = d_CurLocalNodeID;
      MOERTEL::Node n ( d_CurLocalNodeID , pos , 1 , dof , onBoundary , 0 );
      d_CurLocalNodeID++;
      interface.AddNode ( n , side );
      begin++;
    }
  }

  size_t  MoertelOperatorBuilder::buildInformation ( AMP::Mesh::MeshManager::Adapter::shared_ptr mesh
                                                 , short int bid
                                                 , AMP::Mesh::DOFMap::shared_ptr dofs
                                                 , MOERTEL::Interface &interface
                                                 , int side 
                                                 , size_t n_offset 
                                                 , size_t s_offset )
  {
    getNodes ( mesh , mesh->beginOwnedBoundary ( bid ) , mesh->endOwnedBoundary ( bid ) , interface , dofs , side , n_offset );
    int id = 1;
    SideIterator curSide = mesh->beginSideBoundary ( bid );
    while ( curSide != mesh->endSideBoundary ( bid ) )
    {
      int *t = new int [curSide->numNodes()];
      for ( size_t i = 0 ; i != curSide->numNodes() ; i++ )
      {
        t[i] = d_NodeRenumber[curSide->getNodeID(i) + n_offset];
      }
      MOERTEL::Segment_BiLinearQuad s ( id+s_offset , curSide->numNodes () , t , 0 );
      interface.AddSegment ( s , side );
      delete [] t;
      id++;
      curSide++;
    }
    return id;
  }

  void MoertelOperatorBuilder::buildLambdaVector ()
  {
    size_t numLambda = d_M->NumGlobalRows();
    AMP::LinearAlgebra::VectorEngine::BufferPtr  p ( new AMP::LinearAlgebra::VectorEngine::Buffer ( numLambda ) );
    AMP::LinearAlgebra::VectorEngineParameters::shared_ptr  ve_params ( new AMP::LinearAlgebra::EpetraVectorEngineParameters ( d_M->NumMyRows () , d_M->NumGlobalRows () , d_Comm.Comm() ) );

    AMP::LinearAlgebra::VectorEngine::shared_ptr newEngine ( new AMP::LinearAlgebra::EpetraVectorEngine ( ve_params , p ) );
    AMP::LinearAlgebra::ManagedVectorParameters *params = new AMP::LinearAlgebra::ManagedVectorParameters ();
    params->d_Engine = newEngine;
    params->d_CloneEngine = false;
    params->d_Buffer = p;
    d_LambdaVector = AMP::LinearAlgebra::Vector::shared_ptr ( new AMP::LinearAlgebra::ManagedEpetraVector ( AMP::LinearAlgebra::VectorParameters::shared_ptr ( params ) ) );
    d_LambdaVector->setVariable ( AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::Variable ( "lambda" ) ) );
  }

  MoertelOperatorBuilder::MoertelOperatorBuilder ( boost::shared_ptr<MoertelOperatorBuilderParameters> params )
    :  d_Comm ( params->d_Comm.getCommunicator() )
    ,  d_Manager ( d_Comm , MOERTEL::Manager::manager_3D , 0 )
    ,  d_DB ( params->d_DB )
  {
    d_Manager.SetDimension ( MOERTEL::Manager::manager_3D);
    MOERTEL::Interface interface ( 1 , false , d_Comm , 0 );
    d_CurLocalNodeID = 0;
    size_t s_offset = buildInformation ( params->d_MasterVolume , params->d_MasterSurface , params->d_MasterDOFMap , interface , 0 , 0 , 0 );
    s_offset = buildInformation ( params->d_SlaveVolume , params->d_SlaveSurface , params->d_SlaveDOFMap , interface , 1 , params->d_MasterVolume->numTotalNodes() , s_offset );
    interface.SetMortarSide ( 0 );
    interface.SetFunctionTypes ( MOERTEL::Function::func_BiLinearQuad , MOERTEL::Function::func_DualBiLinearQuad );
    if ( !interface.Complete () )
      AMP_ERROR( "Failed to build MOERTEL interface" );
    d_Manager.AddInterface ( interface );
    d_NumNodes = params->d_MasterVolume->numTotalNodes() +
                 params->d_SlaveVolume->numTotalNodes();
    Epetra_Map rowMap ( d_NumNodes , 0 , d_Comm );
    d_Manager.SetProblemMap ( &rowMap );
    Epetra_CrsMatrix prob ( Copy , rowMap , 1 );
    double val =1.; 
    for ( int i = 0 ; i != (int)d_NumNodes ; i++ )
    {
      prob.InsertGlobalValues ( i , 1 , &val , &i );
    }
    prob.FillComplete();  //Finish identity   
    d_Manager.SetInputMatrix ( &prob );

    Teuchos::ParameterList & moertelParams = d_Manager.Default_Parameters ();
    moertelParams.set ( "exact values at gauss points" , true );
    moertelParams.set ( "number gaussian points 2D" , 12 );

    if ( !d_Manager.Mortar_Integrate () )
    {
      AMP_ERROR( "Mortar integration failed" );
    }

    const Epetra_CrsMatrix *t_in = d_Manager.MakeSaddleProblem();
    size_t numSaddleRows = t_in->NumMyRows();
    size_t numMasterNodes = params->d_MasterVolume->numTotalNodes();
    size_t numSlaveNodes = params->d_SlaveVolume->numTotalNodes();
    subMatrix ( t_in , d_M , d_NumNodes , numSaddleRows , 0 , numMasterNodes );
    subMatrix ( t_in , d_D , d_NumNodes , numSaddleRows , numMasterNodes , numMasterNodes + numSlaveNodes );
    subMatrix ( t_in , d_MT , 0 , numMasterNodes , d_NumNodes , numSaddleRows );
    subMatrix ( t_in , d_DT , numMasterNodes , numMasterNodes + numSlaveNodes , d_NumNodes , numSaddleRows );
    buildLambdaVector ();
  }

  MoertelOperatorBuilder::~MoertelOperatorBuilder ()
  {
    delete d_MT;
    delete d_DT;
    delete d_M;
    delete d_D;
  }

  boost::shared_ptr<Operator>    MoertelOperatorBuilder::createMTOperator ( AMP::LinearAlgebra::Variable::shared_ptr invar , AMP::LinearAlgebra::Variable::shared_ptr outvar )
  {
    boost::shared_ptr<EpetraMatrixOperatorParameters>  params ( new EpetraMatrixOperatorParameters ( d_DB ) );
    params->d_Matrix = d_MT;
    boost::shared_ptr<EpetraMatrixOperator> retVal ( new EpetraMatrixOperator ( params ) );
    retVal->setVariables ( invar , outvar );
    return retVal;
  }

  boost::shared_ptr<Operator>    MoertelOperatorBuilder::createDTOperator ( AMP::LinearAlgebra::Variable::shared_ptr invar , AMP::LinearAlgebra::Variable::shared_ptr outvar )
  {
    boost::shared_ptr<EpetraMatrixOperatorParameters>  params ( new EpetraMatrixOperatorParameters ( d_DB ) );
    params->d_Matrix = d_DT;
    boost::shared_ptr<EpetraMatrixOperator> retVal ( new EpetraMatrixOperator ( params ) );
    retVal->setVariables ( invar , outvar );
    return retVal;
  }

  boost::shared_ptr<Operator>    MoertelOperatorBuilder::createMOperator ( AMP::LinearAlgebra::Variable::shared_ptr invar , AMP::LinearAlgebra::Variable::shared_ptr outvar )
  {
    boost::shared_ptr<EpetraMatrixOperatorParameters>  params ( new EpetraMatrixOperatorParameters ( d_DB ) );
    params->d_Matrix = const_cast<Epetra_CrsMatrix *> ( d_M );
    boost::shared_ptr<EpetraMatrixOperator> retVal ( new EpetraMatrixOperator ( params ) );
    retVal->setVariables ( invar , outvar );
    return retVal;
  }

  boost::shared_ptr<Operator>    MoertelOperatorBuilder::createDOperator ( AMP::LinearAlgebra::Variable::shared_ptr invar , AMP::LinearAlgebra::Variable::shared_ptr outvar )
  {
    boost::shared_ptr<EpetraMatrixOperatorParameters>  params ( new EpetraMatrixOperatorParameters ( d_DB ) );
    params->d_Matrix = const_cast<Epetra_CrsMatrix *> ( d_D );
    boost::shared_ptr<EpetraMatrixOperator> retVal ( new EpetraMatrixOperator ( params ) );
    retVal->setVariables ( invar , outvar );
    return retVal;
  }
*/

}
}

