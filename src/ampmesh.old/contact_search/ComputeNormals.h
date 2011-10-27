#ifndef __included_AMP_ComputeNormals_h
#define __included_AMP_ComputeNormals_h

#include "map"

namespace AMP {
namespace Mesh {

  template <typename MESH , typename ITERATOR_FACE>
  class ComputeNormals
  {
    public:
      enum NormalType {Facet , Node};
      class Normal
      {
        public:
          double  x[3];
          Normal () {}
          Normal ( const Normal &rhs )
          {
            for ( int i = 0 ; i != 3 ; i++ )
              x[i] = rhs.x[i];
          }
          Normal &operator = ( const Normal &rhs )
          {
            for ( int i = 0 ; i != 3 ; i++ )
              x[i] = rhs.x[i];
            return *this;
          }
          Normal &operator = ( double d )
          {
            for ( int i = 0 ; i != 3 ; i++ )
              x[i] = d;
            return *this;
          }
          Normal &operator += ( const Normal &rhs )
          {
            for ( int i = 0 ; i != 3 ; i++ )
              x[i] += rhs.x[i];
            return *this;
          }

          Normal &operator -= ( const Normal &rhs )
          {
            for ( int i = 0 ; i != 3 ; i++ )
              x[i] -= rhs.x[i];
            return *this;
          }

          Normal &operator *= ( double d )
          {
            for ( int i = 0 ; i != 3 ; i++ )
              x[i] *= d;
            return *this;
          }

          Normal &operator /= ( double d )
          {
            for ( int i = 0 ; i != 3 ; i++ )
              x[i] /= d;
            return *this;
          }

          double  norm ()
          {
            double denom = 0;
            for ( int i = 0 ; i != 3 ; i++ )
            {
              denom += x[i] * x[i];
            }
            return sqrt ( denom );
          }

          void  normalize ()
          {
            operator /= ( norm() );
          }

          bool operator < ( const Normal &rhs ) const
          {
            for ( size_t i = 0 ; i != 3 ; i++ )
            {
              if ( x[i] < rhs.x[i] ) return true;
              if ( x[i] > rhs.x[i] ) return false;
            }
            return false;
          }
      };

      typedef  std::map<size_t , Normal>       NormalMap;

      static void  compute ( ITERATOR_FACE begin , ITERATOR_FACE end , NormalMap &ans , typename MESH::shared_ptr mesh , NormalType t );
  };


}
}

#include "ComputeNormals.tmpl.h"

#endif

