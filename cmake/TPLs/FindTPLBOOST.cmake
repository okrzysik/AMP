INCLUDE(TribitsTplDeclareLibraries)


TRIBITS_TPL_DECLARE_LIBRARIES( BOOST
    REQUIRED_HEADERS boost/shared_ptr.hpp
)

# Add the definitions
SET( USE_EXT_BOOST 1 )

