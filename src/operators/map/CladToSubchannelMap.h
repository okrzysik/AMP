//#ifndef included_AMP_CladToSubchannelMap
#define included_AMP_CladToSubchannelMap

#include "operators/map/ScalarZAxisMap.h"


namespace AMP {
namespace Operator {


typedef AMP::Operator::Map3to1to3Parameters  CladToSubchannelMapParameters;


/**
 * \class  CladToSubchannelMap
 * \brief  A class used to map temperature from the clad meshes to the subchannel mesh
 * \details  This class maps a scalar quantity (eg. Temperature) from the outer surface
 *    of the clad meshes of a fuel assembly to a subchannel mesh.
 *    mesh1 - clad meshes
 *    mesh2 - subchannel mesh
 */
class CladToSubchannelMap : public AMP::Operator::AsyncMapOperator
{
public:

    /** \brief  Returns true if MapType = "CladToSubchannelMapParameters"
     * \param[in] s  A string extracted from the MapType line in a MeshToMeshMap db
     * \return  True iff s == "CladToSubchannelMapParameters"
     */
    static bool validMapType ( const std::string &s );

    /** \brief  Typedef to identify the parameters class of this operator
     */
    typedef  CladToSubchannelMapParametersParameters   Parameters;

    //!  The base tag used in communication.
    enum { CommTagBase = 25000 };

    /** \brief   Standard constructor
     * \param[in] params  Input parameters
     */
    CladToSubchannelMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> &params );

    //! Destructor
    virtual ~CladToSubchannelMapParameters ();

    /** \brief   Set a frozen vector for results of the apply operation. 
     * \details  Set a frozen vector for results of the apply operation. 
     * \param result    The results vector
     */
    virtual void  setVector ( AMP::LinearAlgebra::Vector::shared_ptr &result );

    /** \brief   Start a communicative apply operation. 
     * \details  Start a communicative apply operation. 
     */
    virtual void applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a = -1.0, const double b = 1.0);

    /** \brief   Finish a communicative apply operation. 
     * \details  Finish a communicative apply operation. 
     */
    virtual void applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &f,
        const AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr &r,
        const double a = -1.0, const double b = 1.0);


protected:

    // Implimented buildReturn routine
    virtual void buildReturn( AMP::LinearAlgebra::Vector::shared_ptr, const AMP::Mesh::Mesh::shared_ptr,
        const AMP::Mesh::MeshIterator&, const std::map<double,double>& );
};


} // Operator namespace
} // AMP namespace


