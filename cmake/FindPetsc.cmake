FUNCTION ( PETSC_GET_VERSION PETSC_INCLUDE_DIR )
    if (EXISTS "${PETSC_INCLUDE_DIR}/petscversion.h")
        FILE (STRINGS "${PETSC_INCLUDE_DIR}/petscversion.h" vstrings REGEX "#define PETSC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
        foreach (line ${vstrings})
            STRING (REGEX REPLACE " +" ";" fields ${line}) # break line into three fields (the first is always "#define")
            LIST (GET fields 1 var)
            LIST (GET fields 2 val)
            SET (${var} ${val} PARENT_SCOPE)
            SET (${var} ${val})         # Also in local scope so we have access below
        endforeach ()
        if (PETSC_VERSION_RELEASE)
            SET (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}p${PETSC_VERSION_PATCH}" PARENT_SCOPE)
        else ()
            # make dev version compare higher than any patch level of a released version
            SET (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.99" PARENT_SCOPE)
        endif ()
    else ()
        MESSAGE (SEND_ERROR "PETSC_DIR can not be used, ${PETSC_DIR}/include/petscversion.h does not exist")
    endif ()
ENDFUNCTION ()


FUNCTION ( PETSC_SET_LIBRARIES  PETSC_LIB_DIRECTORY )
    IF ( "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}" STREQUAL "3.0" )
        FIND_LIBRARY ( PETSC_LIB           NAMES petsc        PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        FIND_LIBRARY ( PETSC_CONTRIB_LIB   NAMES petsccontrib PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        FIND_LIBRARY ( PETSC_PETSCDM_LIB   NAMES petscdm      PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        FIND_LIBRARY ( PETSC_PETSCKSP_LIB  NAMES petscksp     PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        FIND_LIBRARY ( PETSC_PETSCMAT_LIB  NAMES petscmat     PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        FIND_LIBRARY ( PETSC_PETSCSNES_LIB NAMES petscsnes    PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        FIND_LIBRARY ( PETSC_PETSCTS_LIB   NAMES petscts      PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        FIND_LIBRARY ( PETSC_PETSCVEC_LIB  NAMES petscvec     PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        IF ( (NOT PETSC_LIB) OR (NOT PETSC_CONTRIB_LIB) OR (NOT PETSC_PETSCDM_LIB) OR (NOT PETSC_PETSCKSP_LIB) OR 
            (NOT PETSC_PETSCMAT_LIB) OR (NOT PETSC_PETSCSNES_LIB) OR (NOT PETSC_PETSCTS_LIB) OR (NOT PETSC_PETSCVEC_LIB) )
            MESSAGE ( ${PETSC_LIB} )
            MESSAGE ( ${PETSC_CONTRIB_LIB} )
            MESSAGE ( ${PETSC_PETSCDM_LIB} )
            MESSAGE ( ${PETSC_PETSCKSP_LIB} )
            MESSAGE ( ${PETSC_PETSCMAT_LIB} )
            MESSAGE ( ${PETSC_PETSCSNES_LIB} )
            MESSAGE ( ${PETSC_PETSCTS_LIB} )
            MESSAGE ( ${PETSC_PETSCVEC_LIB} )
            MESSAGE ( FATAL_ERROR "PETsc libraries not found in ${PETSC_LIB_DIRECTORY}" )
        ENDIF ()
        IF ( NOT USE_EXT_MPI ) 
            FIND_LIBRARY ( PETSC_MPIUNI_LIB  NAMES mpiuni  PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            IF ( NOT PETSC_MPIUNI_LIB )
                MESSAGE ( ${PETSC_MPIUNI_LIB} )
                MESSAGE ( FATAL_ERROR "PETsc libraries not found in ${PETSC_LIB_DIRECTORY}" )
            ENDIF ()
        ENDIF()
        # Add the libraries in the appropriate order
        SET ( PETSC_LIBS
            ${PETSC_PETSCSNES_LIB}
            ${PETSC_PETSCKSP_LIB}
            ${PETSC_PETSCDM_LIB}
            ${PETSC_CONTRIB_LIB}
            ${PETSC_PETSCMAT_LIB}
            ${PETSC_PETSCVEC_LIB}
            ${PETSC_LIB}
        )
        IF ( NOT USE_EXT_MPI ) 
            SET ( PETSC_LIBS  ${PETSC_LIBS} ${PETSC_MPIUNI_LIB} )
        ENDIF()
        # Set petsc-hypre info
        #FILE ( GLOB HYPRE_FILES "${PETSC_LIB_DIRECTORY}/*HYPRE*" )
        #IF ( HYPRE_FILES AND NOT HYPRE_DIRECTORY )
        #    SET ( HYPRE_DIRECTORY ${PETSC_DIRECTORY}/${PETSC_ARCH} )
        #    CONFIGURE_HYPRE_LIBRARIES ()
        #    SET ( PETSC_LIBS ${PETSC_LIBS} ${HYPRE_LIBS} )
        #ENDIF ()
    ELSEIF( "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}" STREQUAL "3.2" )
        FIND_LIBRARY ( PETSC_LIB           NAMES petsc        PATHS ${PETSC_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        IF ( (NOT PETSC_LIB)  )
            MESSAGE ( ${PETSC_LIB} )
            MESSAGE ( FATAL_ERROR "PETsc libraries not found in ${PETSC_LIB_DIRECTORY}" )
        ENDIF ()
        SET ( PETSC_LIBS ${PETSC_LIB} )
    ELSE()
        MESSAGE ( FATAL_ERROR "AMP is not tested with this version of petsc (${PETSC_VERSION}), use 3.0 or 3.2" )
    ENDIF()
    SET ( PETSC_LIBS ${PETSC_LIBS} PARENT_SCOPE )
ENDFUNCTION ()


