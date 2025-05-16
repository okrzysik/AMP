#ifndef included_AMP_GhostDataHelper
#define included_AMP_GhostDataHelper

#include "AMP/vectors/CommunicationList.h"
#include "AMP/vectors/data/VectorData.h"


namespace AMP::LinearAlgebra {


template<class TYPE = double, class Allocator = AMP::HostAllocator<void>>
class GhostDataHelper : public VectorData
{
public:
    using ScalarAllocator_t =
        typename std::allocator_traits<Allocator>::template rebind_alloc<TYPE>;

    GhostDataHelper();
    GhostDataHelper( std::shared_ptr<CommunicationList> );
    ~GhostDataHelper();

public: // Functions overloaded from VectorData
    bool hasGhosts() const override { return d_ghostSize > 0; }
    std::shared_ptr<CommunicationList> getCommunicationList() const override;
    void setCommunicationList( std::shared_ptr<CommunicationList> comm ) override;
    void aliasGhostBuffer( std::shared_ptr<VectorData> in ) override;
    size_t getGhostSize() const override;
    void fillGhosts( const Scalar & ) override;
    void setNoGhosts() override;
    bool containsGlobalElement( size_t ) const override;
    void setGhostValuesByGlobalID( size_t, const size_t *, const void *, const typeID & ) override;
    void addGhostValuesByGlobalID( size_t, const size_t *, const void *, const typeID & ) override;
    void getGhostValuesByGlobalID( size_t, const size_t *, void *, const typeID & ) const override;
    void
    getGhostAddValuesByGlobalID( size_t, const size_t *, void *, const typeID & ) const override;
    UpdateState getLocalUpdateStatus() const override;
    void setUpdateStatus( UpdateState state ) override;
    void setUpdateStatusPtr( std::shared_ptr<UpdateState> rhs ) override;
    std::shared_ptr<UpdateState> getUpdateStatusPtr() const override;
    void makeConsistent( ScatterType t ) override;
    void dataChanged() override;
    const AMP_MPI &getComm() const override;
    void dumpGhostedData( std::ostream &out, size_t offset ) const override;
    void copyGhostValues( const VectorData &rhs ) override;

    using VectorData::addGhostValuesByGlobalID;
    using VectorData::getGhostAddValuesByGlobalID;
    using VectorData::getGhostValuesByGlobalID;
    using VectorData::makeConsistent;
    using VectorData::setGhostValuesByGlobalID;

public: // Write/read restart data
    void registerChildObjects( AMP::IO::RestartManager *manager ) const override;
    void writeRestart( int64_t fid ) const override;
    GhostDataHelper( int64_t fid, AMP::IO::RestartManager *manager );

protected:
    void scatter_set();
    void scatter_add();
    void deallocateBuffers();
    void allocateBuffers( size_t len );

protected:
    std::shared_ptr<CommunicationList> d_CommList = nullptr;
    std::shared_ptr<UpdateState> d_UpdateState    = nullptr;

    ScalarAllocator_t d_alloc;
    TYPE *d_Ghosts    = nullptr;
    TYPE *d_AddBuffer = nullptr;
    size_t d_ghostSize = 0; //! size/length of ghost and add buffers
};


} // namespace AMP::LinearAlgebra

#endif
