

namespace AMP {
namespace LinearAlgebra {

  inline
  size_t   ManagedDataMap::firstRow () 
  { 
    return d_FirstRow;
  }

  inline
  ManagedDataMap::iterator   ManagedDataMap::begin() 
  { 
    return d_vMapping->begin(); 
  }

  inline
  ManagedDataMap::iterator   ManagedDataMap::end()   
  { 
    return d_vMapping->end(); 
  }

  inline
  ManagedDataMap::const_iterator  ManagedDataMap::begin() const 
  { 
    return d_vMapping->begin(); 
  }

  inline
  ManagedDataMap::const_iterator  ManagedDataMap::end() const   
  { 
    return d_vMapping->end(); 
  }

  inline
  size_t   ManagedDataMap::getLocalSize () const  
  { 
    return d_uiLocalSize; 
  }

  inline
  size_t   ManagedDataMap::getGlobalSize () const 
  { 
    return d_uiGlobalSize; 
  }

  inline
  ManagedDataMap::ManagedDataMap ( size_t local_size , size_t global_size )
    : d_vMapping ( new mapping ( local_size ) ),
      d_uiLocalSize ( local_size ) ,
      d_uiGlobalSize ( global_size ) 
  {
    d_FirstRow = 0;
    d_FirstRowFilled = false;
  }


  inline
  void ManagedDataMap::addMapping ( size_t local_id , size_t global_id )
  {
    if ( local_id >= d_uiLocalSize ) 
      AMP_ERROR( "local_id out of bounds" );

    (*d_vMapping)[local_id] = global_id;
    if ( d_FirstRowFilled )
    {
      d_FirstRow = std::min ( global_id , d_FirstRow );
    }
    else
    {
      d_FirstRow = global_id;
      d_FirstRowFilled = true;
    }
  }

  inline
  ManagedDataMap::ManagedDataMap ( const ManagedDataMap &m )
  {
    d_vMapping = m.d_vMapping;
    d_uiLocalSize = m.d_uiLocalSize;
    d_uiGlobalSize = m.d_uiGlobalSize;
  }

}
}
