
namespace AMP {
namespace LinearAlgebra {

inline ObjectSorterParameters::ObjectSorterParameters() {}

inline ObjectSorter::ObjectSorter() {}

inline ObjectSorter::ObjectSorter( const ObjectSorter & ) {}

inline ObjectSorter::ObjectSorter( Parameters::shared_ptr params )
{
    d_FirstObject = params->d_FirstObject;
}

inline ObjectSorter::~ObjectSorter() {}

inline size_t ObjectSorter::getFirstObject() { return d_FirstObject; }
} // namespace LinearAlgebra
} // namespace AMP
