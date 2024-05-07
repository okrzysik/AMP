#ifndef included_AMP_RandomSubsetVariable
#define included_AMP_RandomSubsetVariable

#include "AMP/vectors/SubsetVariable.h"


namespace AMP::LinearAlgebra {

class RandomSubsetVariable : public SubsetVariable
{
private:
    VectorIndexer::shared_ptr d_Indexer;

public:
    RandomSubsetVariable( const std::string &name, VectorIndexer::shared_ptr i )
        : SubsetVariable( name ), d_Indexer( i )
    {
    }

    virtual VectorIndexer::shared_ptr getIndexer() { return d_Indexer; }
    virtual size_t DOFsPerObject() const { return 0; }
    std::string className() const override { return "RandomSubsetVariable"; }
    void writeRestart( int64_t ) const override;
    RandomSubsetVariable( int64_t );
};
} // namespace AMP::LinearAlgebra

#endif
