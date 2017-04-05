
#ifndef included_AMP_NoteToGeomType::FaceContactOperator
#define included_AMP_NoteToGeomType::FaceContactOperator


#include "utils/shared_ptr.h"
#include <fstream>
#include <matrices/Matrix.h>
#include <operators/contact/ContactOperator.h>
#include <utils/Utilities.h>
#include <vectors/Variable.h>
#include <vectors/Vector.h>

namespace AMP {
namespace Operator {

/**
  An abstract base class for representing a linear operator. This class
  stores the matrix representation of the linear operator. It provides
  an implementation of the apply() function.
  @see Operator
  */
class NodeToGeomType::FaceContactOperator : public ContactOperator
{

public:
    /**
      Constructor. This resets the matrix shared pointer.
      @param [in] params
      */
    NodeToGeomType::FaceContactOperator( const AMP::shared_ptr<ContactOperatorParameters> &params )
        : ContactOperator( params ), d_ContactIsFrictionless( false )
    {
        size_t rank          = d_GlobalComm.getRank();
        std::string fileName = "debug_operator_" + AMP::Utilities::intToString( rank );
        d_fout.open( fileName.c_str(), std::fstream::out );
    }

    /**
      Destructor
      */
    virtual ~NodeToGeomType::FaceContactOperator() { d_fout.close(); }

    /**
     * This function is useful for re-initializing/updating an operator
     * \param params
     *        parameter object containing parameters to change
     */
    void reset( const AMP::shared_ptr<OperatorParameters> &params );

    void addSlaveToMaster( AMP::LinearAlgebra::Vector::shared_ptr u );

    void copyMasterToSlave( AMP::LinearAlgebra::Vector::shared_ptr u );

    void setSlaveToZero( AMP::LinearAlgebra::Vector::shared_ptr u );

    void addShiftToSlave( AMP::LinearAlgebra::Vector::shared_ptr u );

    void initialize();

    size_t updateActiveSet( AMP::LinearAlgebra::Vector::shared_ptr displacementFieldVector,
                            bool skipDisplaceMesh = false );

    size_t updateActiveSetWithALittleHelp( AMP::LinearAlgebra::Vector::shared_ptr helpVector )
    {
        AMP_ASSERT( d_ActiveSet.empty() ); // meant to be use as initial guess only!
        AMP::LinearAlgebra::Vector::shared_ptr masterDisplacement = helpVector->select(
            AMP::LinearAlgebra::VS_Mesh( d_Mesh->Subset( d_MasterMeshID ) ), "dummy" );
        AMP_ASSERT( masterDisplacement->L2Norm() ==
                    0.0 ); // need to change slightly correction otherwise
        d_Mesh->displaceMesh( helpVector );
        size_t const holdOn          = updateActiveSet( helpVector, true );
        size_t const sizeOfActiveSet = d_ActiveSet.size();
        std::vector<size_t> dofIndices;
        double shiftCorrection[3];
        for ( size_t i = 0; i < sizeOfActiveSet; ++i ) {
            d_DOFManager->getDOFs( d_ActiveSet[i], dofIndices );
            AMP_ASSERT( dofIndices.size() == 3 );
            helpVector->getLocalValuesByGlobalID( 3, &( dofIndices[0] ), &( shiftCorrection[0] ) );
            std::transform( &( d_SlaveShift[3 * i] ),
                            &( d_SlaveShift[3 * i] ) + 3,
                            &( shiftCorrection[0] ),
                            &( d_SlaveShift[3 * i] ),
                            std::plus<double>() );
        } // end for i
        helpVector->scale( -1.0 );
        d_Mesh->displaceMesh( helpVector );
        return holdOn;
    }

    void uglyHack( AMP::LinearAlgebra::Vector::shared_ptr temperatureFieldVector,
                   AMP::Discretization::DOFManager::shared_ptr temperatureDOFManager,
                   double thermalExpansionCoefficient,
                   double referenceTemperature )
    {
        d_TemperatureFieldVector      = temperatureFieldVector;
        d_TemperatureDOFManager       = temperatureDOFManager;
        d_ThermalExpansionCoefficient = thermalExpansionCoefficient;
        d_ReferenceTemperature        = referenceTemperature;
    }

    void getSlaveVerticesNormalVectorAndSurfaceTraction(
        std::vector<double> const *&normalVector,
        std::vector<double> const *&surfaceTraction ) const
    {
        normalVector    = &d_SlaveVerticesNormalVectorBeforeUpdate;
        surfaceTraction = &d_SlaveVerticesSurfaceTractionBeforeUpdate;
    }

    void getSlaveVerticesNormalVectorAndShift( std::vector<double> const *&normalVector,
                                               std::vector<double> const *&shift ) const
    {
        normalVector = &d_SlaveVerticesNormalVector;
        shift        = &d_SlaveShift;
    }

    void setContactIsFrictionless( bool isItReally ) { d_ContactIsFrictionless = isItReally; }

protected:
private:
    AMP::LinearAlgebra::Vector::shared_ptr d_TemperatureFieldVector;
    AMP::Discretization::DOFManager::shared_ptr d_TemperatureDOFManager;
    double d_ThermalExpansionCoefficient;
    double d_ReferenceTemperature;

    void getVectorIndicesFromGlobalIDs( const std::vector<AMP::Mesh::MeshElementID> &globalIDs,
                                        std::vector<size_t> &vectorIndices );

    std::vector<int> d_SendCnts;
    std::vector<int> d_SendDisps;
    std::vector<int> d_RecvCnts;
    std::vector<int> d_RecvDisps;
    std::vector<int> d_TransposeSendCnts;
    std::vector<int> d_TransposeSendDisps;
    std::vector<int> d_TransposeRecvCnts;
    std::vector<int> d_TransposeRecvDisps;

    // actually we don't need to store the meshelementids but maybe useful later to check whether
    // the active set has
    // changed
    std::vector<AMP::Mesh::MeshElementID> d_SlaveVerticesGlobalIDs;
    std::vector<AMP::Mesh::MeshElementID> d_MasterVerticesGlobalIDs;
    std::vector<AMP::Mesh::MeshElementID> d_RecvMasterVerticesGlobalIDs;

    //        std::vector<size_t> d_SlaveIndices;
    std::vector<size_t> d_RecvMasterIndices;
    std::vector<size_t> d_MasterVerticesMap;
    std::vector<size_t> d_MasterVerticesInverseMap;
    std::vector<size_t> d_MasterVerticesOwnerRanks;
    std::vector<double> d_MasterShapeFunctionsValues;
    //        std::vector<double> d_SlaveVerticesShift;

    std::vector<AMP::Mesh::MeshElementID> d_MasterVolumesGlobalIDs;
    std::vector<size_t> d_MasterFacesLocalIndices;
    std::vector<double> d_SlaveVerticesNormalVector;

    std::vector<double> d_SlaveVerticesSurfaceTractionBeforeUpdate;
    std::vector<double> d_SlaveVerticesNormalVectorBeforeUpdate;

    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_InputVariable;  /**< Input variable */
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_OutputVariable; /**< Output variable */

    bool d_ContactIsFrictionless;

    std::fstream d_fout;
};

struct ProjectionData {
    AMP::Mesh::MeshElementID d_MasterVolumeGlobalID;
    size_t d_MasterFaceLocalIndex;
    double d_SlaveGeomType::VertexLocalCoordOnMasterFace[2];
};

struct StressStateData {
    double d_SlaveGeomType::VertexNormalVector[3];
    double d_SlaveGeomType::VertexSurfaceTraction[3];
};

struct AnotherDataWithNoName {
    double d_NormalVector[3];
    double d_Displacement[3];
};

struct GeomType::FaceData {
    AMP::Mesh::MeshElementID d_GeomType::FaceVerticesGlobalIDs[4];
};
}
}

#endif
