#include "utils/shared_ptr.h"
#include <vector>
#include <list>
#include <complex>
#include <string>
#include "Database.h"
#include "DatabaseBox.h"
#include "MemoryDatabase.h"
#include "Parser.h"

template class std::list< AMP::shared_ptr<AMP::Database> >;
template class std::vector<std::complex<double> >;
template class std::vector< std::string >;
template class std::vector<unsigned char >;
template class std::vector<AMP::DatabaseBox >;
//template class AMP::DebugAppender::shared_ptr;
//template class AMP::Database::shared_ptr;
//template class AMP::MemoryDatabase::shared_ptr;
//template class boost::detail::sp_counted_impl_p<AMP::MemoryDatabase>;


/*
 * These explicit template instatitions do not work when _GLIBCXX_DEBUG and 
 * GLIBCXX_DEBUG_PEDANTIC are defined.
 */

// template class std::_List_base<AMP::Parser::ParseData, std::allocator<AMP::Parser::ParseData> >;
// template class std::_List_base<AMP::MemoryDatabase::KeyData, std::allocator<AMP::MemoryDatabase::KeyData> >;
// template std::_List_iterator<AMP::MemoryDatabase::KeyData> std::list<AMP::MemoryDatabase::KeyData, std::allocator<AMP::MemoryDatabase::KeyData> >::erase(std::_List_iterator<AMP::MemoryDatabase::KeyData>);
