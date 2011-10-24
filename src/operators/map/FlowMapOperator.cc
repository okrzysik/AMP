#include "FlowMapOperator.h"


namespace AMP {
namespace Operator {


// applyStart
void FlowMapOperator::applyStart ( const AMP::LinearAlgebra::Vector::shared_ptr & , const AMP::LinearAlgebra::Vector::shared_ptr &u ,
        AMP::LinearAlgebra::Vector::shared_ptr & , const double , const double )
{
    d_1DMultiMap.clear();
    if ( d_FirstApply )
    {
      d_FirstApply = false;
    }
    else
    {
      waitForAllRequests ();
    }

    buildMap ( u->subsetVectorForVariable ( getInputVariable() ) , d_TempFlow);
    buildMap ( d_ResultVector , d_TempClad );

    d_SendBuf.resize ( 4 * d_TempFlow.size() );
    std::vector<double>::iterator curBuf = d_SendBuf.begin();
    std::multimap<double,double>::iterator curData = d_TempFlow.begin();
    while ( curData != d_TempFlow.end() )
    {
      *curBuf = curData->first;
      curBuf++;
      *curBuf = curData->second;
      curBuf++;
      curData++;
    }
    curData = d_TempClad.begin();
    while ( curData != d_TempClad.end() )
    {
      *curBuf = curData->first;
      curBuf++;
      *curBuf = curData->second;
      curBuf++;
      curData++;
    }

    std::vector<MPI_Request>::iterator  curReq = beginRequests ();
    for ( size_t i = 0 ; i != d_SendToProc.size() ; i++ )
    {
      if ( d_SendToProc[i] > 0 )
      {
        *curReq = d_MapComm.Isend( getBufferToAvoidDebugVectorCrashing(d_SendBuf), d_SendBuf.size(), i, d_SendTag );
        curReq++;
      }
    }
}


// applyFinish
void  FlowMapOperator::applyFinish ( const AMP::LinearAlgebra::Vector::shared_ptr &f , const AMP::LinearAlgebra::Vector::shared_ptr &u ,
                      AMP::LinearAlgebra::Vector::shared_ptr &r , const double a , const double b )
{
    d_1DMultiMap.clear();
    d_TempClad.clear();
    d_TempFlow.clear();

    for ( size_t i = 0 ; i != d_SendToProc.size() ; i++ )
    {
      if ( d_SendToProc[i] > 0 )
      {
        int inSize = d_MapComm.probe(i,d_RecvTag)/sizeof(double);
        std::vector<double> recvBuf ( inSize );
        int length = recvBuf.size();
        d_MapComm.recv( getBufferToAvoidDebugVectorCrashing( recvBuf ), length, i, false, d_RecvTag );
        int k = inSize/2;
        for ( int j = 0 ; j != inSize/4 ; j++ )
        {
          d_TempFlow.insert ( std::pair<double,double> ( recvBuf[2*j] , recvBuf[2*j+1] ) );
          d_TempClad.insert ( std::pair<double,double> ( recvBuf[2*j+k] , recvBuf[2*j+1+k] ) );
        }
      }
    }
    smear ( 1.e-8 , d_TempFlow );
    smear ( 1.e-8 , d_TempClad );

    std::multimap<double,double>::iterator  T_c_i , T_b_i , T_b_im1;
    T_b_i = T_b_im1 = d_TempFlow.begin();
    T_c_i = d_TempClad.begin();
    T_b_i++;
    T_c_i++;
    addTo1DMap ( T_b_im1->first , 300.0 );
    while ( T_b_i != d_TempFlow.end() )
    {
      AMP_ASSERT ( T_c_i != d_TempClad.end() );

        double Heff, he_z;
        double R_b =0;

        Heff = (0.023*d_K/d_De)*pow(d_Re,0.8)*pow(d_Pr,0.4);
     //       Cp   = getHeatCapacity(T_b_i);
        he_z = T_b_i->first - T_b_im1->first;

        R_b = T_b_i->second - T_b_im1->second - ((4*Heff*( T_c_i->second - T_b_i->second))/(Cp*d_G*d_De))*he_z;
      addTo1DMap ( T_b_i->first , R_b );

      T_b_i++;
      T_b_im1++;
      T_c_i++;
    }


    AMP::LinearAlgebra::Vector::shared_ptr outputVec = r->subsetVectorForVariable ( d_MapVariable );
    buildReturn ( outputVec );

    outputVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

    if(f.get() == NULL) {
        outputVec->scale(a);
    } else {
        AMP::LinearAlgebra::Vector::shared_ptr fInternal = f->subsetVectorForVariable( d_MapVariable );
        if(fInternal.get() == NULL) {
            outputVec->scale(a);
        } else {
            outputVec->axpby(b, a, fInternal);
        }
    }

    d_SendBuf.resize ( 0 );
}


// buildReturn
void FlowMapOperator::buildReturn ( AMP::LinearAlgebra::Vector::shared_ptr p )
{
        AMP::Mesh::DOFMap::shared_ptr d = d_MeshAdapter->getDOFMap ( getInputVariable() );
        BNIterator cur = d_BeginNode;

        if ( d_BeginNode != d_EndNode )
          AMP_ASSERT ( d_1DMultiMap.size() );

        while ( cur != d_EndNode )
        {
          std::multimap<double,double>::iterator lb,ub,be;
          lb = d_1DMultiMap.lower_bound ( cur->z() );
          ub = d_1DMultiMap.upper_bound ( cur->z() );
          be = d_1DMultiMap.begin();
          if ( lb == d_1DMultiMap.end() )
          {
            lb--;
            if ( fabs ( lb->first - cur->z() ) > 1.e-6 )
            {
              cur++;
              continue;
            }
          }

          if ( fabs ( lb->first - cur->z() ) < 1.e-12 )
          {
            double ans = lb->second;
            p->setValueByGlobalID ( d->getGlobalID ( cur->globalID () , 0 ) , ans );
            cur++;
            continue;
          }

          if ( ( lb == d_1DMultiMap.begin() ) || ( ub == d_1DMultiMap.end() ) )
          {
            cur++;
            continue;
          }

          lb--;
          double lo = lb->first;
          double hi = ub->first;
          double pos = cur->z();
          double wt = (pos - lo) / (hi - lo);
          double ans = (1.-wt) * lb->second + wt * ub->second;
          p->setValueByGlobalID ( d->getGlobalID ( cur->globalID() , 0 ) , ans );

          cur++;
        }
}


// buildMap
void FlowMapOperator::buildMap ( const AMP::LinearAlgebra::Vector::shared_ptr p )
{
    AMP_ERROR("This version of buildMap should not be called");
}


// buildMap
void FlowMapOperator::buildMap ( const AMP::LinearAlgebra::Vector::shared_ptr p , std::multimap<double,double> &m )
{
    AMP::Mesh::DOFMap::shared_ptr d = d_MeshAdapter->getDOFMap ( getInputVariable() );
    BNIterator cur = d_BeginNode;
    while ( cur != d_EndNode ) {
        double val = p->getValueByGlobalID ( d->getGlobalID ( cur->globalID() , 0 ) );
        m.insert ( std::pair<double,double> ( cur->z() , val ) );
        cur++;
    }
}


// smear
void FlowMapOperator::smear ( double tolerance )
{
    AMP_ERROR("This version of buildMap should not be called");
}


// smear
void FlowMapOperator::smear ( double tolerance , std::multimap<double,double> &m )
{
    std::multimap<double,double>::iterator  curPtr = m.begin();
    while ( curPtr != m.end() )
    {
      std::multimap<double,double>::iterator  curLookAhead = curPtr;
      while ( curLookAhead != m.end() )
      {
        double cp = curPtr->first;
        double cl = curLookAhead->first;
        if ( fabs ( cp - cl ) > tolerance )
          break;
        curLookAhead++;
      }
      if ( curPtr != curLookAhead )
      {
        double total , count;
        total = count = 0.;
        std::multimap<double,double>::iterator  temp = curPtr;
        while ( temp != curLookAhead )
        {
          total += temp->second;
          count += 1.0;
          temp++;
        }
        curPtr->second = total / count;
        temp = curPtr;
        temp++;
        m.erase ( temp , curLookAhead );
      }
      curPtr++;
    }
}


} // Operators
} // AMP

