#
# Define the list of TPLs, their find module names, and their classification
#
# TPL_NAME:
#
#   The name of the TPL used in the CMake cache variables TPL_ENABLE_${TPL_NAME}
#
# TPL_FINDMOD:
#
#   The name of the find module under that is used to get the names of the
#   TPLs.  If ends in '/' then this gives the directory and the standard module
#   name will be used which is FindTPL${TPL_NAME}.cmake.
#
# TPL_CLASSIFICATION:
#
#   PS: Primary Stable TPL
#
#     Primary Stable TPLs are those TPLs that a developer must have
#     installed on their machine in order to be able to do development.  
#
#   SS: Secondary Stable TPL
#
#     Secondary Stable TPLs are those TPLs that are not required in order to
#     be able to develop and test before checkins but are offically supported.
#     Support for SS TPLs is tested as part of the nightly testing process.
#
#   TS: Tertiary Stable TPL
#
#     Tertiary Stable TPLs are those TPLs that are supported TPLs but can not
#     be included in the set of SS TPLs because they may conflicit with other
#     SS Code.  
#
#   EX: Experimental TPL
#
#     Experimental TPLs are not offically supported.  
#     They represent experimental capabilities.  
#
# The default enable for all TPLs is empty "" reguardless of the category.
# The idea is that the enabling of the TPL will be done by the package and
# other enables that the user has to set.
#

# Get the source directory if it is not set yet
IF ( NOT AMP_SOURCE_DIR )
    STRING(REGEX REPLACE "/PackagesList.cmake" "" AMP_SOURCE_DIR ${CMAKE_CURRENT_LIST_FILE} )
ENDIF()

# Set the TPLs
SET(AMP_TPLS_FINDMODS_CLASSIFICATIONS 
    MPI                 "${AMP_SOURCE_DIR}/cmake/TPLs/"         SS
    BLAS                "${AMP_SOURCE_DIR}/cmake/TPLs/"         PS
    LAPACK              "${AMP_SOURCE_DIR}/cmake/TPLs/"         PS
    BOOST               "${AMP_SOURCE_DIR}/cmake/TPLs/"         PS
#    ZLIB                "${AMP_SOURCE_DIR}/cmake/TPLs/"         PS
    HDF5                "${AMP_SOURCE_DIR}/cmake/TPLs/"         SS
    SILO                "${AMP_SOURCE_DIR}/cmake/TPLs/"         SS
    SUNDIALS            "${AMP_SOURCE_DIR}/cmake/TPLs/"         SS
    PETSC               "${AMP_SOURCE_DIR}/cmake/TPLs/"         SS
    PETSC_AMP           "${AMP_SOURCE_DIR}/cmake/TPLs/"         SS
    LIBMESH             "${AMP_SOURCE_DIR}/cmake/TPLs/"         SS
    DENDRO              "${AMP_SOURCE_DIR}/cmake/TPLs/"         SS
    X11                 "${AMP_SOURCE_DIR}/cmake/TPLs/"         SS
)


