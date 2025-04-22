#ifndef included_AMP_MultiDOFHelper
#define included_AMP_MultiDOFHelper

#include "AMP/discretization/DOF_Manager.h"

#include <vector>

namespace AMP::LinearAlgebra {
class VectorData;
}


namespace AMP::Discretization {

/**
 * \class multiDOFHelper
 * \brief A class to manage mapping multiple DOFs to global indicies
 * \details  This class provides mapping between local and global DOFs for multiple DOF managers
 */
class multiDOFHelper
{
public:
    // Constructor
    multiDOFHelper( const std::vector<std::shared_ptr<DOFManager>> &managers,
                    const AMP::AMP_MPI &comm );

    // Constructor
    multiDOFHelper( const std::vector<AMP::LinearAlgebra::VectorData *> &data,
                    const AMP::AMP_MPI &comm );

    //! Convert the local to global dof
    size_t subToGlobal( int manager, size_t dof ) const;

    //! Convert the local to global dof
    std::vector<size_t> getSubDOF( const int manager, const std::vector<size_t> &globalDOFs ) const;
    //! Convert the global to local dof
    std::pair<size_t, int> globalToSub( size_t dof ) const;

    //! Convert the global to local dof
    std::vector<size_t> getGlobalDOF( const int manager, const std::vector<size_t> &subDOFs ) const;

    // Get the number of local DOFs
    inline size_t numLocal() const { return d_N_local; }

    // Get the number of global DOFs
    inline size_t numGlobal() const { return d_N_global; }

    // Get the first local dof
    inline size_t begin() const { return d_begin; }

    // Get one past the last local dof
    inline size_t end() const { return d_begin + d_N_local; }


public: // Default constructors
    multiDOFHelper()                                    = default;
    multiDOFHelper( multiDOFHelper && )                 = default;
    multiDOFHelper( const multiDOFHelper & )            = default;
    multiDOFHelper &operator=( multiDOFHelper && )      = default;
    multiDOFHelper &operator=( const multiDOFHelper & ) = default;


public: // HDF5 interfaces
    void writeHDF5( size_t ) const;
    multiDOFHelper( size_t );

public:
    // Data used to convert between the local (sub) and global (parent) DOFs
    struct DOFMapStruct {
        // Constructors
        inline DOFMapStruct( size_t sub_start, size_t sub_end, size_t global_start, size_t id )
        {
            data[0] = sub_start;
            data[1] = sub_end;
            data[2] = global_start;
            data[3] = id;
        }
        inline DOFMapStruct()
        {
            data[0] = 0;
            data[1] = 0;
            data[2] = 0;
            data[3] = 0;
        }
        // Convert ids
        inline size_t toGlobal( size_t local ) const { return local - data[0] + data[2]; }
        inline size_t toLocal( size_t global ) const { return global - data[2] + data[0]; }
        inline bool inRangeLocal( size_t local ) const
        {
            return local >= data[0] && local < data[1];
        }
        inline size_t inRangeGlobal( size_t global ) const
        {
            return global >= data[2] && ( global - data[2] ) < ( data[1] - data[0] );
        }
        inline size_t id() const { return data[3]; }
        inline bool operator==( const DOFMapStruct &rhs )
        {
            return data[0] == rhs.data[0] && data[1] == rhs.data[1] && data[2] == rhs.data[2] &&
                   data[3] == rhs.data[3];
        }

    private:
        size_t data[4];
    };

private:
    size_t d_N_local  = 0;
    size_t d_N_global = 0;
    size_t d_begin    = 0;
    std::vector<size_t> d_ids;
    std::vector<DOFMapStruct> d_dofMap;
};


} // namespace AMP::Discretization

#endif
