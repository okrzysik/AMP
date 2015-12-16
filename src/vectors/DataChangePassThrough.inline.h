namespace AMP {
namespace LinearAlgebra {

inline DataChangePassThrough::DataChangePassThrough() {}

inline DataChangePassThrough::~DataChangePassThrough() {}

inline void DataChangePassThrough::dataChanged() { fireDataChange(); }
}
}
