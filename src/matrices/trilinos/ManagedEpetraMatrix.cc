#include "ManagedEpetraMatrix.h"
#include "vectors/trilinos/EpetraVector.h"
#include "vectors/trilinos/EpetraVectorEngine.h"
#include "vectors/trilinos/ManagedEpetraVector.h"
#include <algorithm>
#include "utils/Utilities.h"

#include "utils/AMP_MPI.h"

#include "EpetraExt_MatrixMatrix.h"
#ifdef USE_MPI
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

namespace AMP {
namespace LinearAlgebra {

void ManagedEpetraMatrix::multiply ( shared_ptr other_op , shared_ptr &result )
{
    if ( this->numGlobalColumns() != other_op->numGlobalRows() )
        AMP_ERROR( "Inner matrix dimensions must agree" );
    if ( !other_op->isA<ManagedEpetraMatrix>() )
        AMP_ERROR( "Incompatible matrix types" );

    #ifdef USE_MPI
        MPI_Comm epetraComm = (dynamic_cast<const Epetra_MpiComm *> (&d_epetraMatrix->RowMap().Comm()))->Comm();
    #else
        MPI_Comm epetraComm = AMP_COMM_SELF;
    #endif
    Vector::shared_ptr leftVec = this->getRightVector();
    Vector::shared_ptr rightVec = other_op->getRightVector();
    ManagedEpetraMatrixParameters *memp = new ManagedEpetraMatrixParameters ( 
                 d_epetraMatrix->RowMap().NumMyElements() , 
                 d_epetraMatrix->RowMap().NumGlobalElements() , 
                 d_epetraMatrix->RowMap().MinMyGID() , 
                 AMP_MPI(epetraComm) );
    memp->d_CommListLeft = leftVec->getCommunicationList();
    memp->d_CommListRight = rightVec->getCommunicationList();
    memp->d_DOFManagerLeft = leftVec->getDOFManager();
    memp->d_DOFManagerRight = rightVec->getDOFManager();
    ManagedEpetraMatrix *res = new ManagedEpetraMatrix ( MatrixParameters::shared_ptr ( memp ) );
    EpetraExt::MatrixMatrix::Multiply ( *d_epetraMatrix , false , *(other_op->castTo<ManagedEpetraMatrix> ().d_epetraMatrix) , false , *(res->d_epetraMatrix) , true );
    result = Matrix::shared_ptr ( res );
}


Vector::shared_ptr ManagedEpetraMatrix::getRightVector ()
{
    ManagedEpetraMatrixParameters &memp = d_pParameters->castTo<ManagedEpetraMatrixParameters> ();
    int  localSize = memp.d_CommListRight->numLocalRows();
    int  globalSize = memp.d_CommListRight->getTotalSize();
    EpetraVectorEngineParameters *evep = new EpetraVectorEngineParameters  ( localSize , 
                                                                             globalSize ,
                                                                             memp.getEpetraRowMapPtr () ,
                                                                             memp.getEpetraComm() );
    VectorEngineParameters::shared_ptr p_eng ( evep );
    ManagedVectorParameters *p_params = new ManagedVectorParameters;
    p_params->d_Buffer = VectorEngine::BufferPtr ( new std::vector<double> ( localSize ) );
    p_params->d_Engine = VectorEngine::shared_ptr ( new EpetraVectorEngine ( p_eng , p_params->d_Buffer ) );
    p_params->d_CommList = memp.d_CommListRight;
    p_params->d_DOFManager = memp.d_DOFManagerRight;
    Vector::shared_ptr rtn = Vector::shared_ptr( new ManagedEpetraVector ( VectorParameters::shared_ptr ( p_params ) ) );
    rtn->setVariable( memp.d_VariableRight );
    //rtn->setVariable( Variable::shared_ptr( new Variable("right") ) );
    return rtn;
}


Vector::shared_ptr ManagedEpetraMatrix::getLeftVector ()
{
    ManagedEpetraMatrixParameters &memp = d_pParameters->castTo<ManagedEpetraMatrixParameters> ();
    int  localSize = memp.d_CommListLeft->numLocalRows();
    int  globalSize = memp.d_CommListLeft->getTotalSize();
    EpetraVectorEngineParameters *evep = new EpetraVectorEngineParameters  ( localSize , 
                                                                             globalSize ,
                                                                             memp.getEpetraColMapPtr () ,
                                                                             memp.getEpetraComm() );
    VectorEngineParameters::shared_ptr p_eng ( evep );
    ManagedVectorParameters *p_params = new ManagedVectorParameters;
    p_params->d_Buffer = VectorEngine::BufferPtr ( new std::vector<double> ( localSize ) );
    p_params->d_Engine = VectorEngine::shared_ptr ( new EpetraVectorEngine ( p_eng , p_params->d_Buffer ) );
    p_params->d_CommList = memp.d_CommListLeft;
    p_params->d_DOFManager = memp.d_DOFManagerLeft;
    Vector::shared_ptr rtn = Vector::shared_ptr( new ManagedEpetraVector ( VectorParameters::shared_ptr ( p_params ) ) );
    rtn->setVariable( memp.d_VariableLeft );
    //rtn->setVariable( Variable::shared_ptr( new Variable("left") ) );
    return rtn;
}


Discretization::DOFManager::shared_ptr ManagedEpetraMatrix::getRightDOFManager ()
{
    ManagedEpetraMatrixParameters &memp = d_pParameters->castTo<ManagedEpetraMatrixParameters> ();
    return memp.d_DOFManagerRight;
}


Discretization::DOFManager::shared_ptr ManagedEpetraMatrix::getLeftDOFManager ()
{
    ManagedEpetraMatrixParameters &memp = d_pParameters->castTo<ManagedEpetraMatrixParameters> ();
    return memp.d_DOFManagerLeft;
}


void ManagedEpetraMatrix::setOtherData ()
{
    AMP_MPI myComm(d_pParameters->castTo<ManagedEpetraMatrixParameters>().getEpetraComm());
    int ndxLen = d_OtherData.size();
    int totNdxLen = myComm.sumReduce(ndxLen);
    if ( totNdxLen == 0 ) {
        return;
    }
    int  dataLen = 0;
    std::map<int,std::map<int , double> >::iterator cur_row = d_OtherData.begin();
    while ( cur_row != d_OtherData.end() ) {
        dataLen += cur_row->second.size();
        cur_row++;
    }
    int *rows = new int [ dataLen+1 ];   //Add one to have the new work
    int *cols = new int [ dataLen+1 ];
    double *data = new double [ dataLen+1 ];
    int cur_ndx = 0;
    int cur_ptr = 0;
    cur_row = d_OtherData.begin();
    while ( cur_row != d_OtherData.end() )
    {
      cur_ndx++;
      std::map<int,double>::iterator cur_elem = cur_row->second.begin();
      while ( cur_elem != cur_row->second.end() )
      {
        rows[cur_ptr] = cur_row->first;
        cols[cur_ptr] = cur_elem->first;
        data[cur_ptr] = cur_elem->second;
        cur_ptr++;
        cur_elem++;
      }
      cur_row++;
    }

    int totDataLen = myComm.sumReduce(dataLen);

    int *aggregateRows = new int [ totDataLen ];
    int *aggregateCols = new int [ totDataLen ];
    double *aggregateData = new double [ totDataLen ];

    myComm.allGather( rows, dataLen, aggregateRows );
    myComm.allGather( cols, dataLen, aggregateCols );
    myComm.allGather( data, dataLen, aggregateData );

    int MyFirstRow = d_pParameters->castTo<ManagedEpetraMatrixParameters>().d_DOFManagerLeft->beginDOF();
    int MyEndRow = d_pParameters->castTo<ManagedEpetraMatrixParameters>().d_DOFManagerLeft->endDOF();
    for ( int i = 0 ; i != totDataLen ; i++ )
    {
      if ( ( aggregateRows[i] >= MyFirstRow ) && ( aggregateRows[i] < MyEndRow ) )
      {
        setValueByGlobalID ( aggregateRows[i] , aggregateCols[i] , aggregateData[i] );
      }
    }

    d_OtherData.clear();
    delete [] rows;
    delete [] cols;
    delete [] data;
    delete [] aggregateRows;
    delete [] aggregateCols;
    delete [] aggregateData;
}


ManagedEpetraMatrixParameters::ManagedEpetraMatrixParameters ( int local_size , int global_size , int first_dof , AMP_MPI comm )
        : MatrixParameters () ,
          d_comm ( comm ) ,
          d_vEntriesPerRow ( local_size ) ,
          d_ColGlobal ( global_size ) ,
          d_ColBase ( -1 ) ,
          d_RowBase ( first_dof )
{
}


ManagedEpetraMatrixParameters::ManagedEpetraMatrixParameters ( int local_size , int global_size , int first_dof , int col_global , int col_base , AMP_MPI comm )
        : MatrixParameters () ,
          d_comm ( comm ) ,
          d_vEntriesPerRow ( local_size ) ,
          d_ColGlobal ( col_global ) ,
          d_ColBase ( col_base ) ,
          d_RowBase ( first_dof )
{
}


Epetra_Map  &ManagedEpetraMatrixParameters::getEpetraRowMap ()
{
    #ifdef USE_MPI
        Epetra_MpiComm  comm = d_comm.getCommunicator();
    #else
        Epetra_SerialComm  comm;
    #endif
    AMP_ASSERT(d_DOFManagerLeft.get()!=NULL);
    AMP_ASSERT(d_DOFManagerRight.get()!=NULL);
    size_t N_row_local = d_DOFManagerLeft->numLocalDOF();
    size_t N_row_global = d_DOFManagerLeft->numGlobalDOF();
    if ( d_eRowMap.get() == 0 )
        d_eRowMap = boost::shared_ptr<Epetra_Map>( new Epetra_Map ( N_row_global, N_row_local, d_RowBase, comm ) );
    return *d_eRowMap;
}


Epetra_Map  *ManagedEpetraMatrixParameters::getEpetraColMap ()
{
    #ifdef USE_MPI
        Epetra_MpiComm  comm = d_comm.getCommunicator();
    #else
        Epetra_SerialComm  comm;
    #endif
    if ( d_ColBase < 0 ) return 0;
    std::vector<int> cols;
    cols.reserve ( d_sColumns.size() );
    for ( std::set<int>::iterator curCol = d_sColumns.begin(); curCol != d_sColumns.end() ; curCol++ )
      cols.push_back ( *curCol );
    if ( d_eColMap.get() == 0 )
      d_eColMap = boost::shared_ptr<Epetra_Map> ( new Epetra_Map ( -1 , cols.size() , &*(cols.begin()) , d_ColBase , comm ) );
    return d_eColMap.get();
}


void ManagedEpetraMatrixParameters::addColumns ( int a , int *b )
{
    for ( int i = 0 ; i != a ; i++ )
      d_sColumns.insert ( b[i] );
}


ManagedEpetraMatrix::ManagedEpetraMatrix ( MatrixParameters::shared_ptr params )
        : EpetraMatrix ( params->castTo<ManagedEpetraMatrixParameters>().getEpetraRowMap() , 
                         params->castTo<ManagedEpetraMatrixParameters>().getEpetraColMap()  ,
                         params->castTo<ManagedEpetraMatrixParameters>().entryList() ) ,
          ManagedMatrix ( params ) ,
          d_pParameters ( boost::static_pointer_cast<ManagedEpetraMatrixParameters>(params) )
{
}


ManagedEpetraMatrix::ManagedEpetraMatrix ( const ManagedEpetraMatrix &rhs )
        : Matrix () ,
          EpetraMatrix ( rhs.d_pParameters->castTo<ManagedEpetraMatrixParameters>().getEpetraRowMap() , 
                           rhs.d_pParameters->castTo<ManagedEpetraMatrixParameters>().getEpetraColMap()  ,
                           rhs.d_pParameters->castTo<ManagedEpetraMatrixParameters>().entryList() ) ,
          ManagedMatrix ( rhs.d_pParameters ) ,
          d_pParameters ( rhs.d_pParameters )
{
    ManagedEpetraMatrixParameters  &params = d_pParameters->castTo<ManagedEpetraMatrixParameters> ();
    for ( size_t i = params.d_DOFManagerLeft->beginDOF(); i!=params.d_DOFManagerLeft->endDOF(); i++ ) {
      std::vector<unsigned int>  cols;
      std::vector<double>        vals;
      rhs.getRowByGlobalID ( (int)i , cols , vals );
      for ( size_t j = 0 ; j != cols.size() ; j++ )
      {
        vals[j] = 0;
      }
      if ( cols.size() )
        createValuesByGlobalID ( 1 , (int)cols.size() , (int *)&i , (int *)&(cols[0]) , &(vals[0]) );
    }
    d_RangeMap = rhs.d_RangeMap;
    d_DomainMap = rhs.d_DomainMap;
    makeConsistent();
}


/*
Matrix::shared_ptr  ManagedEpetraMatrix::cloneMatrix () const
{
    return Vector::shared_ptr ( new ManagedEpetraMatrix ( d_pParameters ) );
}
*/


void ManagedEpetraMatrix::mult ( const Vector::shared_ptr &in , Vector::shared_ptr &out )
{
    AMP_ASSERT ( in->getGlobalSize() == numGlobalRows() );
    AMP_ASSERT ( out->getGlobalSize() == numGlobalColumns() );
    Vector::shared_ptr  in_view = EpetraVector::view ( in );
    Vector::shared_ptr  out_view = EpetraVector::view ( out );
    VerifyEpetraReturn ( d_epetraMatrix->Multiply ( false , in_view->castTo<EpetraVector>().getEpetra_Vector() , out_view->castTo<EpetraVector>().getEpetra_Vector() ) , "mult" );
}


void ManagedEpetraMatrix::multTranspose ( const Vector::shared_ptr &in , Vector::shared_ptr &out )
{
    Vector::shared_ptr  in_view = EpetraVector::view ( in );
    Vector::shared_ptr  out_view = EpetraVector::view ( out );
    VerifyEpetraReturn ( d_epetraMatrix->Multiply ( true , in_view->castTo<EpetraVector>().getEpetra_Vector() , out_view->castTo<EpetraVector>().getEpetra_Vector() ) , "mult" );
}


Vector::shared_ptr  ManagedEpetraMatrix::extractDiagonal ( Vector::shared_ptr v )
{
    Vector::shared_ptr  retVal;
    if ( v )
    {
      retVal = EpetraVector::view ( v );
    }
    else
    {
      retVal = getRightVector();
    }
    VerifyEpetraReturn ( d_epetraMatrix->ExtractDiagonalCopy ( retVal->castTo<EpetraVector>().getEpetra_Vector() ) , 
                         "extractDiagonal" );
    return retVal;
}
  

void  ManagedEpetraMatrix::scale ( double alpha ) 
{
    VerifyEpetraReturn ( d_epetraMatrix->Scale ( alpha ) , "scale" );
}


void  ManagedEpetraMatrix::axpy ( double alpha , const Matrix &rhs )
{
    int    *a , *b;
    double *values1;
    double *values2;
    d_epetraMatrix->ExtractCrsDataPointers ( a , b , values1 );
    rhs.castTo<ManagedEpetraMatrix>().d_epetraMatrix->ExtractCrsDataPointers ( a , b , values2 );
    for ( int i = 0 ; i != d_epetraMatrix->NumMyNonzeros() ; i++ )
      values1[i] += alpha*values2[i];
}


void  ManagedEpetraMatrix::addValuesByGlobalID ( int  num_rows , int num_cols , int *rows , int *cols , double *values )
{
    for ( int i = 0 ; i != num_rows ; i++ )
      VerifyEpetraReturn ( d_epetraMatrix->SumIntoGlobalValues ( rows[i] , num_cols , values+num_cols*i , cols ) , "addValuesByGlobalId" );
}


void  ManagedEpetraMatrix::createValuesByGlobalID ( int  num_rows , int num_cols , int *rows , int *cols , double *values )
{
    for ( int i = 0 ; i != num_rows ; i++ )
      VerifyEpetraReturn (d_epetraMatrix->InsertGlobalValues ( rows[i] , num_cols , values+num_cols*i , cols ) , "setValuesByGlobalID" );
}


void  ManagedEpetraMatrix::setValuesByGlobalID ( int  num_rows , int num_cols , int *rows , int *cols , double *values )
{
    ManagedEpetraMatrixParameters &EpetraParamters = d_pParameters->castTo<ManagedEpetraMatrixParameters>();

    int MyFirstRow = EpetraParamters.d_DOFManagerLeft->beginDOF();
    int MyEndRow = EpetraParamters.d_DOFManagerLeft->endDOF();
    for ( int i = 0 ; i != num_rows ; i++ )
    {
      VerifyEpetraReturn (d_epetraMatrix->ReplaceGlobalValues ( rows[i] , num_cols , values+num_cols*i , cols ) , "setValuesByGlobalID" );
      for ( int j = 0 ; j != num_cols ; j++ )
      {
        if ( (values[num_cols*i + j] < MyFirstRow ) ||
             (values[num_cols*i + j] >= MyEndRow ) )
        {
          d_OtherData[rows[i]][cols[j]] = values[num_cols*i + j];
        }
      }
    }
}


void  ManagedEpetraMatrix::setScalar ( double ans )
{
    VerifyEpetraReturn ( d_epetraMatrix->PutScalar ( ans ) , "setScalar" );
}


void  ManagedEpetraMatrix::makeConsistent ()
{
    Epetra_FECrsMatrix  *mat = dynamic_cast<Epetra_FECrsMatrix *> ( d_epetraMatrix );
    if ( mat )
    {
      VerifyEpetraReturn ( mat->GlobalAssemble ( false ) , "makeParallelConsistent" );
      fillComplete ();
    }
    setOtherData ();
}


void ManagedEpetraMatrix::setDiagonal ( const Vector::shared_ptr &in )
{
    VerifyEpetraReturn ( d_epetraMatrix->ReplaceDiagonalValues ( EpetraVector::view ( in )->castTo<EpetraVector>().getEpetra_Vector() ) , "setDiagonal" );
}


void ManagedEpetraMatrix::getRowByGlobalID ( int row , std::vector<unsigned int> &cols , std::vector<double> &values ) const
{   
    int firstRow = d_pParameters->d_DOFManagerLeft->beginDOF();
    int numRows = d_pParameters->d_DOFManagerLeft->endDOF();
    AMP_ASSERT ( row >= firstRow );
    AMP_ASSERT ( row < firstRow + numRows );

    int localRow = row - firstRow;
    int numCols = d_pParameters->entriesInRow ( localRow );
    cols.resize ( numCols );
    values.resize ( numCols );

    if ( numCols )
      VerifyEpetraReturn ( d_epetraMatrix->ExtractGlobalRowCopy ( row , numCols , numCols , &(values[0]) , (int *)&(cols[0]) ) , "getRowByGlobalID" );
}


}
}

