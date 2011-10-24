#ifndef included_AMP_MeshUtils
#define included_AMP_MeshUtils

#include "utils/AMP_MPI.h"


namespace AMP { 
namespace Mesh {


template <typename OPERATOR>
typename OPERATOR::result_type  MeshAccumulate ( typename OPERATOR::iterator begin ,
                                                 typename OPERATOR::iterator end ,
                                                 OPERATOR op = OPERATOR () )
{
    typename OPERATOR::result_type  accumulator = op.zero ( begin );
    while ( begin != end )
    {
      op ( begin , accumulator );
      begin++;
    }
    return accumulator;
}


template <typename OPERATOR>
void   MeshApply ( typename OPERATOR::iterator begin , typename OPERATOR::iterator end ,
                   OPERATOR op = OPERATOR () )
{
    while ( begin != end )
    {
      op ( begin );
      begin++;
    }
}


template <typename OPERATOR , typename CONSTRUCTABLE>
void   MeshConstruct ( typename OPERATOR::iterator begin , typename OPERATOR::iterator end ,
                       std::vector<CONSTRUCTABLE> &vec , OPERATOR op = OPERATOR() )
{
    vec.clear();
    while ( begin != end )
    {
      vec.push_back ( op ( begin ) );
      begin++;
    }
}


struct simple_point
{
    double x , y , z;
};


template <typename T>
struct min_max_struct
{
    T  min;
    T  max;
};


template <typename ITERATOR>
struct  min_max_point
{
    typedef min_max_struct<simple_point>  result_type;
    typedef ITERATOR                      iterator;
    
    result_type  zero ( ITERATOR  begin )
    {
      result_type retval;
      retval.min.x = retval.max.x = begin->x();
      retval.min.y = retval.max.y = begin->y();
      retval.min.z = retval.max.z = begin->z();
      return retval;
    }

    void  operator ()( ITERATOR  cur , result_type &accumulator )
    {
      accumulator.min.x = std::min ( accumulator.min.x , cur->x() );
      accumulator.min.y = std::min ( accumulator.min.y , cur->y() );
      accumulator.min.z = std::min ( accumulator.min.z , cur->z() );
      accumulator.max.x = std::max ( accumulator.max.x , cur->x() );
      accumulator.max.y = std::max ( accumulator.max.y , cur->y() );
      accumulator.max.z = std::max ( accumulator.max.z , cur->z() );
    }
};


template <typename ITERATOR>
struct  min_max_radius
{
    typedef  min_max_struct<double>     result_type;
    typedef ITERATOR                    iterator;
    simple_point                        middle;

    double      compute_radius ( ITERATOR cur )
    {
      return sqrt ( cur->x()*cur->x() + cur->y()*cur->y() );
    }

    result_type zero ( ITERATOR  begin )
    {
      result_type retval;
      retval.min = retval.max = compute_radius ( begin );
      return retval;
    }

    void  operator ()( ITERATOR cur , result_type &accumulator )
    {
      accumulator.min = std::min ( accumulator.min , compute_radius ( cur ) );
      accumulator.max = std::max ( accumulator.max , compute_radius ( cur ) );
    }
};


template <typename ADAPTER>
min_max_struct<simple_point>  computeExtremeCoordinates ( typename ADAPTER::shared_ptr &adapter )
{
    min_max_struct<simple_point>  local , retval;
    local = MeshAccumulate<min_max_point<typename ADAPTER::NodeIterator> > ( adapter->beginNode() , adapter->endNode() );
    adapter->getComm().minReduce((double*)&local.min,(double*)&retval.min,3);
    adapter->getComm().maxReduce((double*)&local.max,(double*)&retval.max,3);
    return retval;
}


template <typename ADAPTER>
min_max_struct<double>  computeExtremeRadii ( typename ADAPTER::shared_ptr adapter )
{
    min_max_struct<double>  local , retval;
    local = MeshAccumulate<min_max_radius<typename ADAPTER::NodeIterator> > ( adapter->beginNode() , adapter->endNode() );
    adapter->getComm().minReduce((double*)&local.min,(double*)&retval.min,1);
    adapter->getComm().maxReduce((double*)&local.max,(double*)&retval.max,1);
    return retval;
}


} // Mesh namespace
} // AMP namespace


#endif
