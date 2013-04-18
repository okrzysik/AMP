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

SET(AMP_TPLS_FINDMODS_CLASSIFICATIONS 
  MPI               "cmake/TPLs/"                           SS
  BLAS              "cmake/TPLs/"                           PS
  LAPACK            "cmake/TPLs/"                           PS
  BOOST             "cmake/TPLs/"                           PS
  HDF5              "cmake/TPLs/"                           SS
  SILO              "cmake/TPLs/"                           SS
  Zlib              "cmake/TPLs/"                           SS
  PETSC             "cmake/TPLs/"                           SS
  LIBMESH           "cmake/TPLs/"                           SS
  TRILINOS          "cmake/TPLS/"                           SS
)
SET(SCALE_TPLS_FINDMODS_CLASSIFICATIONS ${Scale_TPLS_FINDMODS_CLASSIFICATIONS} )

