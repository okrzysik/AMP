#include "utils/AMP_MPI.h"

namespace AMP { 
namespace Mesh {

  template <typename ADAPTER>
  NodalRowMapTmplParameters<ADAPTER>::NodalRowMapTmplParameters ( typename ADAPTER::ElementIterator begin,
                                                                  typename ADAPTER::ElementIterator end )
    : d_beginElement ( begin )
    , d_endElement ( end )
  {
  }

//  template <typename ADAPTER>
//  NodalRowMapTmpl<ADAPTER>::NodalRowMapTmpl () 
//  {
//  }

  template <typename ADAPTER>
  NodalRowMapTmpl<ADAPTER>::~NodalRowMapTmpl () 
  {
  }

  template <typename ADAPTER>
  unsigned int *NodalRowMapTmpl<ADAPTER>::getRowIndices ()
  {
    return &(d_vRowNdxs[0]);
  }

  template <typename ADAPTER>
  unsigned int *NodalRowMapTmpl<ADAPTER>::getColumns ( int row )
  {
    return &(d_vColumns[d_vRowNdxs[row]]);
  }

  template <typename ADAPTER>
  unsigned int NodalRowMapTmpl<ADAPTER>::getNNZ ( int row )
  {
    return d_vRowNdxs[row+1] - d_vRowNdxs[row];
  }

//  template <typename ADAPTER>
//  unsigned int NodalRowMapTmpl<ADAPTER>::numLocalRows () const
//  {
//    return d_vRowNdxs.size() - 1;
//  }

  template <typename ADAPTER>
  void NodalRowMapTmpl<ADAPTER>::distributeRemoteGraph ( Graph &local , Graph &remote )
  {
    size_t  totLen = 1;  // Once for total number of nodes;
    NodeIterator  cur_rem_node = remote.begin();
    while ( cur_rem_node != remote.end() )
    {
      totLen += 2 + cur_rem_node->second.size();
      // Once for the node id
      // Once for total number of edges;
      // Once for each adjacent node;
      cur_rem_node++;
    }

    int  *graphSerial = new int [ totLen ];
    size_t  curSpot = 0;
    graphSerial[curSpot] = remote.size();
    curSpot++;
    cur_rem_node = remote.begin();
    while ( cur_rem_node != remote.end() )
    {
      graphSerial[curSpot] = cur_rem_node->first;
      curSpot++;
      graphSerial[curSpot] = cur_rem_node->second.size();
      curSpot++;
      NeighborIterator  cur_neighbor = cur_rem_node->second.begin();
      while ( cur_neighbor != cur_rem_node->second.end() )
      {
        graphSerial[curSpot] = *cur_neighbor;
        curSpot++;
        cur_neighbor++;
      }
      cur_rem_node++;
    }

    AMP_MPI myComm = AMP_MPI( d_params->d_mesh->getComm() );

    int  aggregateLen = myComm.sumReduce(totLen);
    int  *aggregateGraphs = new int [ aggregateLen ];

    myComm.allGather ( graphSerial , totLen , aggregateGraphs );
    curSpot = 0;
    while ( (int)curSpot < aggregateLen )
    {
      size_t  numNodes = aggregateGraphs[curSpot];
      curSpot++;
      for ( size_t i = 0 ; i != numNodes ; i++ )
      {
        NodeIterator myNode = local.find ( aggregateGraphs[curSpot] );
        curSpot++;
        size_t  numNeighbors = aggregateGraphs[curSpot];
        curSpot++;
        for ( size_t j = 0 ; j != numNeighbors ; j++ )
        {
          if ( myNode != local.end() )
          {
            myNode->second.insert ( aggregateGraphs[curSpot] );
          }
          curSpot++;
        }
      }
    }

    delete [] graphSerial;
    delete [] aggregateGraphs;
  }


template <typename ADAPTER>
void NodalRowMapTmpl<ADAPTER>::buildSparsityPattern ( ElementIterator begin , ElementIterator end ,
                                             DOFMap &dofMap , std::vector<unsigned int>  &DOFs )
{
    Graph graph , remote_graph;
    std::set<unsigned int>  dofs;

    // Construct the DOF graphs (this is a slow loop)
    /* The bulk of the time is spent inserting into graph and remote_graph
     * Because these are both std::map, the loop scales as O(N*log(N)*K*log(N)), where 
     * N is the number of DOFs (per processor), and K is the average number of 
     * DOFs connected to each DOF.  Additionally, there are O(N*K) memory allocations.
     */
    size_t dofMap_begin = dofMap.beginDOF();
    size_t dofMap_end = dofMap.endDOF();
    std::vector<unsigned int>  dof_list(100,0);     // Initialize vector storage for the dof list
    for (ElementIterator cur_elem=begin; cur_elem!=end; cur_elem++) {
        // Get the DOFs for the current element
        dof_list.clear();
        dofMap.getDOFs( *cur_elem, dof_list );
        // Sort the DOFs
        Utilities::quicksort(dof_list);
        std::vector<unsigned int>::iterator first = dof_list.begin();
        std::vector<unsigned int>::iterator last  = dof_list.end();
        // For each DOF, add the list of DOFs to the map
        for (size_t i=0; i<dof_list.size(); i++) {
            unsigned int cur_dof = dof_list[i];
            dofs.insert( cur_dof );
            if ( ( cur_dof >= dofMap_begin ) && ( cur_dof < dofMap_end ) ) {
                graph[cur_dof].insert(first,last);
            } else {
                remote_graph[cur_dof].insert(first,last);
            }
        }
    }

    distributeRemoteGraph ( graph , remote_graph );
    int  num_dofs = 0;
    int  num_nz = 0;
    NodeIterator  cur_dof_node = graph.begin();
    while ( cur_dof_node != graph.end() )
    {
      num_dofs++;
      num_nz += cur_dof_node->second.size();
      cur_dof_node++;
    }

    d_vRowNdxs.resize ( num_dofs + 1 );
    d_vColumns.resize ( num_nz );

    int  col_vec_ndx = 0;
    int  row_ndx = 0;
    cur_dof_node = graph.begin();
    while ( cur_dof_node != graph.end() )
    {
      d_vRowNdxs[row_ndx] = col_vec_ndx;
      NeighborIterator cur_neighbor = cur_dof_node->second.begin();
      while ( cur_neighbor != cur_dof_node->second.end() )
      {
        d_vColumns[col_vec_ndx] = *cur_neighbor;
        cur_neighbor++;
        col_vec_ndx++;
      }
      cur_dof_node++;
      row_ndx++;
    }
    d_vRowNdxs[row_ndx] = col_vec_ndx;

    DOFs.resize ( dofs.size() );
    std::set<unsigned int>::iterator  cur_dof = dofs.begin();
    int dof_off = 0;
    while ( cur_dof != dofs.end() )
    {
      DOFs[dof_off] = *cur_dof;
      cur_dof++;
      dof_off++;
    }
  }

  template <typename ADAPTER>
  NodalRowMapTmpl<ADAPTER>::NodalRowMapTmpl ( typename Parameters::shared_ptr  params )
       : CommunicationList ( params )
       , d_params ( params )
  {
  }

  template <typename ADAPTER>
  void NodalRowMapTmpl<ADAPTER>::finalizeList ( Castable &dofMap )
  {
    std::vector<unsigned int> DOFs;
    buildSparsityPattern ( d_params->d_beginElement , d_params->d_endElement , dofMap.castTo<DOFMap>() , DOFs );

    // Build communication map
    AMP_MPI comm = d_params->d_mesh->getComm();
    int commSize = comm.getSize();
    int commRank = comm.getRank();

    std::vector<unsigned int> partitionInfo ( commSize );
    unsigned int t = dofMap.castTo<DOFMap>().endElement();
    comm.allGather(t,&(partitionInfo[0]));

    buildCommunicationArrays ( DOFs , partitionInfo , commRank );
  }


}
}

