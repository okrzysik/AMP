#ifndef included_AMP_RandomAccessIndexer
#define included_AMP_RandomAccessIndexer

#include <map>

#include "AMP/vectors/VectorIndexer.h"


namespace AMP::LinearAlgebra {

class RandomAccessIndexer : public VectorIndexer
{
private:
    // Temporary fix for const correctness.
    mutable std::map<size_t, size_t> d_SuperToSub;
    mutable std::map<size_t, size_t> d_SubToSuper;

public:
    RandomAccessIndexer() {}
    void addID( size_t i ) { d_SuperToSub[i] = 0; }
    void finalize()
    {
        std::map<size_t, size_t>::iterator c = d_SuperToSub.begin();
        size_t lid                           = 0;
        while ( c != d_SuperToSub.end() ) {
            c->second         = lid;
            d_SubToSuper[lid] = c->first;
            lid++;
            c++;
        }
    }

    virtual bool isInSub( size_t i ) const { return d_SuperToSub.find( i ) != d_SuperToSub.end(); }
    virtual size_t getIncrement( size_t i ) const
    {
        std::map<size_t, size_t>::iterator c = d_SuperToSub.find( i );
        AMP_ASSERT( c != d_SuperToSub.end() );
        c++;
        if ( c != d_SuperToSub.end() )
            return ( c->first - i );
        return 2000000000; // For a SubsetVector, the iterator returns end when the
                           // internal iterator iterates past the end;
    }

    virtual size_t getSubID( size_t i ) const { return d_SuperToSub[i]; }
    virtual size_t getSuperID( size_t i ) const { return d_SubToSuper[i]; }
    virtual size_t getNumLocalElements( Vector::shared_ptr v ) const
    {
        size_t start                         = v->getCommunicationList()->getStartGID();
        size_t end                           = v->getCommunicationList()->numLocalRows() + start;
        size_t i                             = 0;
        std::map<size_t, size_t>::iterator c = d_SuperToSub.begin();
        while ( c != d_SuperToSub.end() ) {
            if ( c->first >= start )
                break;
            i++;
            c++;
        }
        std::map<size_t, size_t>::reverse_iterator d = d_SuperToSub.rbegin();
        while ( d != d_SuperToSub.rend() ) {
            if ( d->first < end )
                break;
            i++;
            d++;
        }
        return d_SuperToSub.size() - i;
    }
};
} // namespace AMP::LinearAlgebra

#endif
