#include "utils/shared_ptr.h"
#include "Database.h"
#include "DatabaseBox.h"
#include "MemoryDatabase.h"
#include "Parser.h"
#include <complex>
#include <list>
#include <string>
#include <vector>

template class std::list<AMP::shared_ptr<AMP::Database>>;
template class std::vector<std::complex<double>>;
template class std::vector<std::string>;
template class std::vector<unsigned char>;
template class std::vector<AMP::DatabaseBox>;
