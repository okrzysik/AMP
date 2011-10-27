#ifndef included_AMP_MeshIterators
#define included_AMP_MeshIterators

#include <iterator>

#include "utils/Castable.h"

namespace AMP { 
namespace Mesh {

  template <typename ITERATOR , typename OBJECT_WRAPPER>
  class  MeshIteratorWrapper
  {
    private:
      ITERATOR          d_Iterator;
      mutable OBJECT_WRAPPER    d_Object;

    public:
      // This makes the class a legitimate iterator in the eyes of C++
      typedef  size_t                              difference_type;
      typedef  OBJECT_WRAPPER                      value_type;
      typedef  OBJECT_WRAPPER *                    pointer;
      typedef  OBJECT_WRAPPER &                    reference;
      typedef  std::bidirectional_iterator_tag     iterator_category;

      typedef ITERATOR         Iterator;
      typedef OBJECT_WRAPPER   Object;

      const Iterator       &getIterator () const { return d_Iterator; }
      const Object         &getObject   () const { return d_Object; }

      MeshIteratorWrapper () {}
      MeshIteratorWrapper ( const ITERATOR &rhs ) : d_Iterator ( rhs )
                                                  , d_Object () {}
      template <typename IN_ITERATOR>
      MeshIteratorWrapper ( const MeshIteratorWrapper<IN_ITERATOR,OBJECT_WRAPPER> &rhs )
                                                             : d_Iterator ( rhs.getIterator() )
                                                             , d_Object ( rhs.getObject() ) {}

      MeshIteratorWrapper<ITERATOR,OBJECT_WRAPPER>  & operator++ ()      { ++d_Iterator; return *this;}
      MeshIteratorWrapper<ITERATOR,OBJECT_WRAPPER>    operator++ ( int ) { ++d_Iterator; return *this;}
      MeshIteratorWrapper<ITERATOR,OBJECT_WRAPPER>  & operator-- ()      { --d_Iterator; return *this;}
      MeshIteratorWrapper<ITERATOR,OBJECT_WRAPPER>    operator-- ( int ) { --d_Iterator; return *this;}


      template <typename IN_ITERATOR>
      bool operator == ( const MeshIteratorWrapper<IN_ITERATOR,OBJECT_WRAPPER> &rhs ) const 
                                              { 
                                                Iterator i = rhs.getIterator();
                                                return d_Iterator == i; 
                                              }

      template <typename IN_ITERATOR>
      bool operator != ( const MeshIteratorWrapper<IN_ITERATOR,OBJECT_WRAPPER> &rhs ) const 
                                              { 
                                                Iterator i = rhs.getIterator();
                                                return d_Iterator != i; 
                                              }

      OBJECT_WRAPPER   operator *  () const { return d_Object = *d_Iterator; }
      OBJECT_WRAPPER  &operator *  ()       { return d_Object = *d_Iterator; }
      OBJECT_WRAPPER  *operator -> ()       { return &(d_Object = *d_Iterator); }
      const OBJECT_WRAPPER  *operator -> () const      { return &(d_Object = *d_Iterator); }

      ITERATOR  getBaseIterator () { return d_Iterator; }

      template <typename IN_ITERATOR>
      bool operator < ( const MeshIteratorWrapper<IN_ITERATOR,OBJECT_WRAPPER> &rhs ) const
      {
        const Object t1 = *(getIterator());
        const Object t2 = *(rhs.getIterator());
        return t1.globalID() < t2.globalID();
      //  ITERATOR myIter = getIterator();
      //  ITERATOR rhsIter = getIterator();
      //  return myIter < rhsIter;
      }
  };

  template <typename ITERATOR , typename OBJECT_WRAPPER>
  class  OwnedMeshIteratorWrapper
  {
    private:
      ITERATOR          d_Iterator;
      ITERATOR          d_End;
      OBJECT_WRAPPER    d_Object;

    public:
      typedef ITERATOR         Iterator;
      typedef OBJECT_WRAPPER   Object;

      const Iterator       &getIterator () const { return d_Iterator; }
      const Object         &getObject   () const { return d_Object; }

      OwnedMeshIteratorWrapper () {}
      OwnedMeshIteratorWrapper ( const ITERATOR &rhs , const ITERATOR &end ) 
                                                  : d_Iterator ( rhs )
                                                  , d_End ( end )
                                                  , d_Object () 
      {
        if ( rhs == end )
          return;
        while ( !(*this)->isOwned() )
        {
          ++d_Iterator; 
          if ( d_Iterator == d_End )
            break;
        }
      }

      template <typename IN_ITERATOR>
      OwnedMeshIteratorWrapper ( const OwnedMeshIteratorWrapper<IN_ITERATOR,OBJECT_WRAPPER> &rhs )
                                                             : d_Iterator ( rhs.getIterator() )
                                                             , d_End ( rhs.d_End )
                                                             , d_Object ( rhs.getObject() ) 
      {
      }

      OwnedMeshIteratorWrapper<ITERATOR,OBJECT_WRAPPER>  & operator++ ()      
      { 
        ++d_Iterator;
        if ( d_Iterator != d_End )
        {
          while ( !(*this)->isOwned() )
          {
            ++d_Iterator; 
            if ( d_Iterator == d_End )
              break;
          }
        }
        return *this;
      }

      OwnedMeshIteratorWrapper<ITERATOR,OBJECT_WRAPPER>    operator++ ( int ) 
      { 
        ++d_Iterator;
        if ( d_Iterator != d_End )
        {
          while ( !(*this)->isOwned() )
          {
            ++d_Iterator; 
            if ( d_Iterator == d_End )
              break;
          }
        }
        return *this;
      }


      template <typename IN_ITERATOR>
      bool operator == ( const OwnedMeshIteratorWrapper<IN_ITERATOR,OBJECT_WRAPPER> &rhs ) const 
                                              { 
                                                Iterator i = rhs.getIterator();
                                                return d_Iterator == i; 
                                              }

      template <typename IN_ITERATOR>
      bool operator != ( const OwnedMeshIteratorWrapper<IN_ITERATOR,OBJECT_WRAPPER> &rhs ) const 
                                              { 
                                                Iterator i = rhs.getIterator();
                                                return d_Iterator != i; 
                                              }

      OBJECT_WRAPPER   operator *  () { return d_Object = *d_Iterator; }
      OBJECT_WRAPPER  *operator -> () { return &(d_Object = *d_Iterator); }

      ITERATOR  getBaseIterator () { return d_Iterator; }

      template <typename IN_ITERATOR>
      bool operator < ( const MeshIteratorWrapper<IN_ITERATOR,OBJECT_WRAPPER> &rhs ) const
      {
        return getObject().globalID() < rhs->globalID();
      }
  };


}
}
#endif



