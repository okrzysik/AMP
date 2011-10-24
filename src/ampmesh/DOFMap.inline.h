namespace AMP { 
namespace Mesh {

  inline
  size_t  DOFMap::beginDOF () 
  { 
    return firstElement(); 
  }

  inline
  size_t  DOFMap::endDOF () 
  { 
    return endElement(); 
  }

}
}

