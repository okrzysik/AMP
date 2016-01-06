#ifndef included_Parser_INLINE
#define included_Parser_INLINE

#include "Parser.h"

namespace AMP {


inline Parser *Parser::getParser() { return ( s_default_parser ); }

inline AMP::shared_ptr<Database> &Parser::getScope() { return ( d_scope_stack.front() ); }

inline void Parser::enterScope( const std::string &name )
{
    d_scope_stack.push_front( d_scope_stack.front()->putDatabase( name ) );
}

inline void Parser::leaveScope() { d_scope_stack.pop_front(); }

inline int Parser::getNumberErrors() const { return ( d_errors ); }

inline int Parser::getNumberWarnings() const { return ( d_warnings ); }
}

#endif
