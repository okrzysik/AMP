INCLUDE(TribitsTplDeclareLibraries)
INCLUDE( ${AMP_SOURCE_DIR}/cmake/macros.cmake )


TRIBITS_TPL_DECLARE_LIBRARIES( X11
    REQUIRED_LIBS_NAMES SM ICE X11
)


MESSAGE ( "Using X11" )
MESSAGE ( "   "  ${X11_LIBS} )


# Add the tribits flags
SET( TPL_ENABLE_X11 ON )


# Add the definitions
SET( USE_EXT_X11 1 )
MESSAGE ( "Using X11" )
MESSAGE ( "   "  ${TPL_X11_LIBRARIES} )

