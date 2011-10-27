
#include "DOFMap.h"
#include "MeshElement.h"


namespace AMP { 
namespace Mesh {

DOFMap::DOFMap ( AMP::LinearAlgebra::VectorEntryMap<true>::Parameters::shared_ptr rhs )
    : AMP::LinearAlgebra::VectorEntryMap<true> ( rhs )
{
    getCommunicationList()->finalize ( *this );
}


void DOFMap::getDOFs ( const MeshObject &obj , std::vector <unsigned int> &ans , unsigned int i ) const
{
    const MeshElement &elem = obj.castTo<MeshElement>();
    size_t numNodes = elem.numNodes ();
    size_t numDOFsPerNode = d_Variable->DOFsPerObject();

    if ( i == static_cast<unsigned int> ( -1 ) )
    {
      ans.resize ( numDOFsPerNode * numNodes );
      for ( size_t j = 0 ; j != numNodes ; j++ )
      {
        for ( size_t k = 0 ; k != numDOFsPerNode ; k++ )
        {
          ans[j*numDOFsPerNode + k] = getGlobalID ( elem.getNodeID ( j ) , k );
        }
      }
    }
    else
    {
      ans.resize ( numNodes );
      for ( size_t k = 0 ; k != numNodes ; k++ )
      {
        ans[k] = getGlobalID ( elem.getNodeID ( k ) , i );
      }
    }
}


void DOFMap::getDOFs ( const MeshObject &obj , std::vector <unsigned int> &ans , std::vector <unsigned int> i ) const
{
    if ( i.size() == 0 ) {
        size_t DOFsPerObject = d_Variable->DOFsPerObject();
        i.resize ( DOFsPerObject );
        for ( size_t k=0 ; k!=DOFsPerObject; k++ )
            i[k] = k;
    }
    ans.resize ( i.size() );
    size_t globalID = obj.globalID();
    for ( size_t k=0; k!=i.size(); k++ )
        ans[k] = getGlobalID( globalID , i[k] );
}


}
}

