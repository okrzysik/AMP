#ifndef included_AMP_VectorSelector_h
#define included_AMP_VectorSelector_h

#include "vectors/Vector.h"

#ifdef USE_AMP_MESH
    #include "ampmesh/Mesh.h"
    #include "ampmesh/MeshIterator.h"
#endif


namespace AMP {
namespace LinearAlgebra {


/** \brief A class used by Vector::select and Vector::selectInto to create
  * vectors with particular data.
  * \details  VectorSelector is designed to perform two types of selection:
  * gross and fine.  The isSelected method will determine if a Vector should
  * be considered for fine selection.  The subset method can be used to create 
  * a fine subset of the Vector.  These methods are
  * meant to be used solely in Vector::selectInto ().  Subclasses of Vector
  * are encouraged to call Vector::selectInto () rather than use these methods.
  */
class VectorSelector 
{
public:
    /** \brief Virtual destructor */
    virtual ~VectorSelector ();

    /** \brief Returns true if Vector grossly matches a selection condition 
      * \param[in]  vec  The Vector to match
      * \details Base class defaults to accepting all vectors.
      */
    virtual  bool  isSelected ( Vector::const_shared_ptr vec ) const;

    /** \brief Subset the given vector
      * \param[in]  vec  The Vector to subset
      * \details Base class defaults to returning all data in the vector
      */
    virtual  Vector::shared_ptr  subset ( Vector::shared_ptr p ) const;
};


/** \brief  Select a subvector based on the name of the variable
  * \details  This can be used in the Vector::select or Vector::selectInto interface:
  * \code
       // Create a vector of all data in the simulation
       Vector::shared_ptr   data = meshManager->createVector ( complexVariable );

       // Extract all data called "Temperature"
       Vector::shared_ptr   temperature = data->select ( VS_ByVariableName ( "Temperature" ) , "Temperature" );

       // Add displacement data to a vector results
       data->select ( VS_ByBariableName ( "Displacement" ) , results );
       \endcode
  */
class VS_ByVariableName : public VectorSelector
{
protected:
    std::string    d_VecName;
public:
    /** \brief Constructor
      * \param[in] n  The name of the variable to subset on
      */
    VS_ByVariableName ( std::string  n );

    virtual bool   isSelected ( Vector::const_shared_ptr v ) const;
};


/** \brief  Select a subvector based on the type of the variable
  * \details  This can be used in the Vector::select or Vector::selectInto interface:
  * \code
       // Create a vector of all data in the simulation
       Vector::shared_ptr   data = meshManager->createVector ( complexVariable );

       // Extract all nodal scalar data
       Vector::shared_ptr   reals = data->select ( VS_ByVariableType<NodalScalarVariable> () , "Real Fields" );

       // Add 3-Vector data to a vector results
       data->select ( VS_ByBariableType <Nodal3VectorVariable> () , "Positions" );
       \endcode
*/
template <typename T>
class VS_ByVariableType : public VectorSelector
{
public:
    virtual bool   isSelected ( Vector::const_shared_ptr v ) const;
};


/** \brief  Create a subset based on a stride in the vector
  * \details  This will pick every \f$b\f$th element starting at \f$a\f$ in an vector
  */
class VS_Stride : public VectorSelector
{
public:
    /** \brief Constructor
      * \param[in]  n  The name of the new variable
      * \param[in]  a  The offset to stride
      * \param[in]  b  The length to stride
      */
    VS_Stride ( const std::string &n , size_t a , size_t b );

    virtual  Vector::shared_ptr  subset ( Vector::shared_ptr p ) const;

protected:
    //  Offset to start striding on
    size_t  d_Offset; 
    //  The stride to use
    size_t  d_Stride;
    //  The name of this subset
    std::string  d_Name;
};


#ifdef USE_AMP_MESH
/** \brief  Create a subset based on a mesh 
  * \details  This will select the portion of a vector that is on the given mesh
  */
class VS_Mesh : public VectorSelector
{
public:
    /** \brief Constructor
      * \param[in]  name    The name of the new variable
      * \param[in]  mesh    The desired mesh
      */
    VS_Mesh ( const std::string &name, AMP::Mesh::Mesh::shared_ptr mesh );

    virtual  Vector::shared_ptr  subset ( Vector::shared_ptr p ) const;

protected:
    //  The name of this subset
    std::string  d_Name;
    //  Mesh
    AMP::Mesh::Mesh::shared_ptr  d_mesh;
};


/** \brief  Create a subset based on a mesh iterator
  * \details  This will select the portion of a vector that is on the given mesh iterator
  */
class VS_MeshIterator : public VectorSelector
{
public:
    /** \brief Constructor
      * \param[in]  name        The name of the new variable
      * \param[in]  iterator    The mesh iterator to use
      */
    VS_MeshIterator ( const std::string &name, const AMP::Mesh::MeshIterator &iterator );

    virtual  Vector::shared_ptr  subset ( Vector::shared_ptr p ) const;

protected:
    //  The name of this subset
    std::string  d_Name;
    //  MeshIterator
    const AMP::Mesh::MeshIterator  d_iterator;
};
#endif


}
}

#include "VectorSelector.tmpl.h"

#endif
