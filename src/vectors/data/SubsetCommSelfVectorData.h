#ifndef included_AMP_SubsetCommSelfVectorData
#define included_AMP_SubsetCommSelfVectorData

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/subsetCommSelfDOFManager.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/data/SubsetVectorData.h"
#include "AMP/vectors/data/VectorData.h"

#include <vector>

namespace AMP::LinearAlgebra {


/** \class SubsetCommSelfVectorData
  * \brief This vector is a subset of an AMP Vector
  * \details  Given an AMP Vector, this class will create a view of a subset of the
        vector for AMP_COMM_SELF.
  */
class SubsetCommSelfVectorData : public VectorData
{

public:
    std::string VectorDataName() const override
    {
        return "SubsetCommSelfVectorData of " + d_parentData->VectorDataName();
    }
    size_t numberOfDataBlocks() const override;
    size_t sizeOfDataBlock( size_t i ) const override;

    void addValuesByLocalID( size_t, const size_t *, const void *, const typeID & ) override;
    void setValuesByLocalID( size_t, const size_t *, const void *, const typeID & ) override;
    void getValuesByLocalID( size_t, const size_t *, void *, const typeID & ) const override;
    void putRawData( const void *in, const typeID &id ) override;
    void getRawData( void *out, const typeID &id ) const override;
    typeID getType( size_t block ) const override { return d_parentData->getType( block ); }
    uint64_t getDataID() const override { return d_parentData->getDataID(); }
    size_t sizeofDataBlockType( size_t ) const override { return sizeof( double ); }
    void swapData( VectorData & ) override;
    std::shared_ptr<VectorData> cloneData( const std::string &name = "" ) const override;
    bool hasContiguousData() const override { return numberOfDataBlocks() > 1 ? false : true; }
    SubsetCommSelfVectorData() {}
    explicit SubsetCommSelfVectorData( std::shared_ptr<SubsetVectorParameters> params );

private:
    void *getRawDataBlockAsVoid( size_t i ) override;
    const void *getRawDataBlockAsVoid( size_t i ) const override;

    // Internal data
    std::shared_ptr<VectorData> d_parentData;          // VectorData for the subset
    size_t d_parentLocalStartID;                       // Offset for the parent
    std::vector<size_t> d_SubsetLocalIDToViewGlobalID; // The list of global ID in the parent vector
    std::shared_ptr<AMP::Discretization::subsetCommSelfDOFManager>
        d_DOFManager; // this will contain the subset DOF
};


} // namespace AMP::LinearAlgebra


#endif
