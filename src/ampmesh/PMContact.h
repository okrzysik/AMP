#ifndef included_AMP_PMContact_h
#define included_AMP_PMContact_h

#include <string>

namespace AMP { 
namespace Mesh {

  template <typename MANAGER>
  struct PMContact
  {
    typedef          MANAGER                   MeshManager;
    typedef typename MANAGER::shared_ptr       MeshManagerPtr;
    typedef typename MANAGER::Adapter          MeshAdapter;
    typedef typename MeshAdapter::shared_ptr   MeshAdapterPtr;
    
//    static void  shiftUp ( MeshManagerPtr manager , std::vector<std::string> &prefix );
    static void  restrictInnerRadius ( MeshManagerPtr manager , double r , const std::string &prefix );
    static void  restrictOuterRadius ( MeshManagerPtr manager , double r , const std::string &prefix );
    static double computeShiftUp ( MeshManagerPtr manager , std::vector<std::string> &prefix );
  };

}
}

#include "PMContact.tmpl.h"
#endif
