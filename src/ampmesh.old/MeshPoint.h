#ifndef included_AMP_MeshPoint_h
#define included_AMP_MeshPoint_h

namespace AMP { 
namespace Mesh {

  class MeshPoint 
  {
    public:
      virtual ~MeshPoint () {}

      virtual double   x () const = 0;
      virtual double   y () const = 0;
      virtual double   z () const = 0;

      virtual double  &x () = 0;
      virtual double  &y () = 0;
      virtual double  &z () = 0;

      virtual double  &operator () ( size_t i ) = 0;
      virtual double   operator () ( size_t i ) const = 0;
  };

}
}

#endif
