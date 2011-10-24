#ifndef included_AMP_NodalRowMap_h
#define included_AMP_NodalRowMap_h

#include "vectors/CommunicationList.h"
#include "DOFMap.h"

namespace AMP { 
namespace Mesh {

  template <typename ADAPTER>
  class NodalRowMapTmplParameters : public AMP::LinearAlgebra::CommunicationListParameters
  {
    public:
      typedef boost::shared_ptr<NodalRowMapTmplParameters <ADAPTER> >   shared_ptr;

      typename ADAPTER::ElementIterator    d_beginElement;
      typename ADAPTER::ElementIterator    d_endElement;
               ADAPTER                    *d_mesh;

      NodalRowMapTmplParameters ( typename ADAPTER::ElementIterator , typename ADAPTER::ElementIterator );
  };

  template <typename ADAPTER>
  class NodalRowMapTmpl : public AMP::LinearAlgebra::CommunicationList
  {
    public:
      typedef  NodalRowMapTmplParameters<ADAPTER>     Parameters;

    private:
      NodalRowMapTmpl ();

    protected:
      int     d_FirstDOF;
      int     d_DOFsPerNode;
      int     d_EndDof;


      typedef typename ADAPTER::ElementIterator   ElementIterator;

      std::vector<unsigned int>           d_vColumns;
      std::vector<unsigned int>           d_vRowNdxs;
      unsigned int                        d_iBeginDof;

      typename Parameters::shared_ptr     d_params;

      void  distributeRemoteGraph ( Graph & , Graph & );
      void  buildSparsityPattern ( ElementIterator begin , ElementIterator end , DOFMap & , std::vector<unsigned int> & );

      virtual  void           finalizeList ( Castable &map );

    public:
      NodalRowMapTmpl ( typename Parameters::shared_ptr ptr );
      virtual ~NodalRowMapTmpl ();

      virtual  unsigned int  *getRowIndices ();
      virtual  unsigned int  *getColumns ( int row );
      virtual  unsigned int   getNNZ ( int row );
//      virtual  unsigned int   numLocalRows () const;
  };


}
}

#include "NodalRowMap.tmpl.h"

#endif
