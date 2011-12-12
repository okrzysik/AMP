#ifndef included_AMP_VectorSubsetter
#define included_AMP_VectorSubsetter

#include "SubsetVector.h"
#include "boost/shared_ptr.hpp"
#include "StridedVariable.h"
#include "RandomAccessIndexer.h"
#include "RandomSubsetVariable.h"

namespace AMP {
namespace LinearAlgebra {

  /** \class VectorSubsetter
    * \brief  A class used to subset a vector
    * \details  Given an AMP vector, a subsetter will create a particular subset for you
    *  The base class returns the entire vector when subset is called
    */
  class VectorSubsetter
  {
    public:
      /** \brief Convenience typedef
        */
      typedef boost::shared_ptr <VectorSubsetter>   shared_ptr;

      /** \brief Destructor
        */
      virtual ~VectorSubsetter ();
      virtual Vector::shared_ptr  subset ( Vector::shared_ptr p );
  };


  /** \class VectorStriderSubsetter
    * \brief  Create subsets based on a StridedVariable
    * \see StridedVariable
    */
  class VectorStriderSubsetter : public VectorSubsetter
  {
    protected:
      /** \brief Offset to start
        */
      size_t   d_Offset;

      /** \brief Stride length
        */
      size_t   d_Stride;

      /** \brief Name of new vector
        */
      std::string  d_Name;

    public:

      /** \brief Constructor
        * \param[in] n  name of new vector
        * \param[in] a  offset to start striding
        * \param[in] b  length of stride
        */
      VectorStriderSubsetter ( const std::string &n , size_t a , size_t b );

      /** \brief Destructor
        */
      virtual ~VectorStriderSubsetter ();

      virtual Vector::shared_ptr  subset ( Vector::shared_ptr p );
  };


  class VectorRandomAccessSubsetter : public VectorSubsetter
  {
    protected:
      std::string                d_Name;
      VectorIndexer::shared_ptr  d_RAIndexer;

    public:
      VectorRandomAccessSubsetter ( const std::string &n );

      void   addID ( size_t i ) { d_RAIndexer->castTo<RandomAccessIndexer> ().addID ( i ); }
      virtual ~VectorRandomAccessSubsetter () {}
      virtual Vector::shared_ptr   subset ( Vector::shared_ptr p );
  };

}
}


#endif
