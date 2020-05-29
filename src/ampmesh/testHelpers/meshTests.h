#ifndef included_AMP_MeshTests
#define included_AMP_MeshTests

#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/UnitTest.h"

#ifdef USE_AMP_VECTORS
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/testHelpers/VectorTests.h"
#endif


namespace AMP {
namespace Mesh {


/**
 * \class meshTests
 * \brief A helper class to store/run tests for a mesh
 */
class meshTests
{
public:
    /**
     * \brief Run all mesh based tests
     * \details  This test runs all the mesh-based tests
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void MeshTestLoop( AMP::UnitTest *ut, std::shared_ptr<AMP::Mesh::Mesh> mesh );

    /**
     * \brief Run all mesh geometry based tests
     * \details  This test runs all the geometry-based tests
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void MeshGeometryTestLoop( AMP::UnitTest *ut, std::shared_ptr<AMP::Mesh::Mesh> mesh );


    /**
     * \brief Run all vector based tests
     * \details  This test runs all the vector-based tests
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     * \param[in] fast          Speed up testing by eliminating some of the tests
     */
    static void MeshVectorTestLoop( AMP::UnitTest *ut,
                                    std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                    bool fast = false );


    /**
     * \brief Run all matrix based tests
     * \details  This test runs all the matrix-based tests
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     * \param[in] fast          Speed up testing by eliminating some of the tests
     */
    static void MeshMatrixTestLoop( AMP::UnitTest *ut,
                                    std::shared_ptr<AMP::Mesh::Mesh> mesh,
                                    bool fast = false );

public: // Basic tests
    /**
     * \brief Check basic id info
     * \details  This tests checks some trivial ids
     * \param[in,out] ut        Unit test class to report the results
     */
    static void testID( AMP::UnitTest *ut );


public: // Mesh based tests
    /**
     * \brief Checks the mesh iterators
     * \details  This test performs a series of simple tests on the basic iterators within a mesh
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void MeshIteratorTest( AMP::UnitTest *ut, std::shared_ptr<AMP::Mesh::Mesh> mesh );


    /**
     * \brief Checks the mesh operators
     * \details  This test performs a series of simple tests on the operator== and operator!=
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void MeshIteratorOperationTest( AMP::UnitTest *ut,
                                           std::shared_ptr<AMP::Mesh::Mesh> mesh );

    /**
     * \brief Checks the mesh set operations
     * \details  This test performs a series of simple tests on the union, intersection, and
     * compliment operations.
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void MeshIteratorSetOPTest( AMP::UnitTest *ut, std::shared_ptr<AMP::Mesh::Mesh> mesh );

    /**
     * \brief Checks the number of elements in a mesh
     * \details  This test performs a series of simple tests on numLocalElements, numGlobalElements,
     * and numGhostElements
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void MeshCountTest( AMP::UnitTest *ut, std::shared_ptr<AMP::Mesh::Mesh> mesh );


    /**
     * \brief Run basic tests
     * \details  This test performs some very basic mesh tests
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void MeshBasicTest( AMP::UnitTest *ut, std::shared_ptr<AMP::Mesh::Mesh> mesh );


    /**
     * \brief Check ghost elements
     * \details  This tests checks that all ghost elements are owned by "owner processor"
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void VerifyGhostIsOwned( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


    /**
     * \brief Check boundary ids
     * \details  This tests loops over all boundary ids, checking the iterators
     * \param[out] ut           Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void VerifyBoundaryIDNodeIterator( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );

    /**
     * \brief Check boundary
     * \details  This tests checks the boundary iterators
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void VerifyBoundaryIterator( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


    /**
     * \brief Check block ids
     * \details  This tests loops over the blocks, checking the iterators
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void testBlockIDs( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


    /**
     * \brief Check neighbors
     * \details  Test if we correctly identify the node neighbors
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void getNodeNeighbors( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


    /**
     * \brief Check displacement by a scalar
     * \details  Test if we correctly displace the mesh
     * \param[out] ut           Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void DisplaceMeshScalar( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


    /**
     * \brief Check displacement by a vector
     * \details  Test if we correctly displace the mesh
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void DisplaceMeshVector( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


    /**
     * \brief Check parents
     * \details  Test getting parent elements for each mesh element
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void getParents( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


    /**
     * \brief Check clone
     * \details  Test cloning a mesh
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void cloneMesh( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


    /**
     * \brief Check mesh performance
     * \details Test the performance of some common mesh operations
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void MeshPerformance( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


public: // Old tests
    static void VerifyNodeElemMapIteratorTest( AMP::UnitTest *ut,
                                               AMP::Mesh::Mesh::shared_ptr mesh );
    static void VerifyBoundaryIteratorTest( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );
    static void VerifyElementForNode( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );

public: // Geometry based tests
    /**
     * \brief Basic geometry class tests
     * \details  This runs the geometry only tests
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh contaning geometry to test
     */
    static void TestBasicGeometry( AMP::UnitTest *ut, AMP::Mesh::Mesh::const_shared_ptr mesh );

    /**
     * \brief Checks Geometry::inside
     * \details  This test checks if all points in the mesh are inside the geometry
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void TestInside( AMP::UnitTest *ut, AMP::Mesh::Mesh::const_shared_ptr mesh );

    /**
     * \brief Checks Geometry::inside
     * \details  This test checks the physical-logical-physical transformation
     * \param[in,out] ut        Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    static void TestPhysicalLogical( AMP::UnitTest *ut, AMP::Mesh::Mesh::const_shared_ptr mesh );


public: // Vector based tests
#ifdef USE_AMP_VECTORS
        //! Factory to create a vector from a mesh
    class MeshVectorFactory : public AMP::LinearAlgebra::VectorFactory
    {
    public:
        MeshVectorFactory( AMP::Mesh::Mesh::shared_ptr mesh,
                           AMP::Mesh::GeomType type,
                           int gcw,
                           int dofs_per_node,
                           bool split )
            : d_dofManager( AMP::Discretization::simpleDOFManager::create(
                  mesh, type, gcw, dofs_per_node, split ) ),
              SPLIT( split )
        {
        }
        //! Get the Variable
        virtual AMP::LinearAlgebra::Variable::shared_ptr getVariable() const override
        {
            return std::make_shared<AMP::LinearAlgebra::Variable>( "test vector" );
        }
        //! Get the Vector
        virtual AMP::LinearAlgebra::Vector::shared_ptr getVector() const override
        {
            return AMP::LinearAlgebra::createVector( d_dofManager, getVariable(), SPLIT );
        }
        //! Get the DOFManager
        virtual AMP::Discretization::DOFManager::shared_ptr getDOFMap() const override
        {
            return d_dofManager;
        }
        //! Get the name
        virtual std::string name() const override { return "MeshVectorFactory"; }

    private:
        MeshVectorFactory();
        AMP::Discretization::DOFManager::shared_ptr d_dofManager;
        bool SPLIT;
    };
#endif

    /**
     * \brief Simple nodal vector tests
     * \details Run a series of simple tests on a nodal vector
     * \param[out] ut           Unit test class to report the results
     * \param[in] mesh          Mesh to test
     * \param[in] DOFs          DOF Manager to use
     * \param[in] gcw           Ghost cell width to use
     */
    template<int DOF_PER_NODE, bool SPLIT>
    static void simpleNodalVectorTests( AMP::UnitTest *ut,
                                        AMP::Mesh::Mesh::shared_ptr mesh,
                                        AMP::Discretization::DOFManager::shared_ptr DOFs,
                                        int gcw );


    /**
     * \brief VerifyGetVectorTest
     * \details VerifyGetVectorTest
     * \param[out] ut           Unit test class to report the results
     * \param[in] mesh          Mesh to test
     */
    template<int DOF_PER_NODE, bool SPLIT>
    static void VerifyGetVectorTest( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


public: // Matrix based tests
    template<int DOF_PER_NODE, bool SPLIT>
    static void VerifyGetMatrixTrivialTest( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );

    template<int DOF_PER_NODE, bool SPLIT>
    static void GhostWriteTest( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh );


private: // Private data
         /**
          * \brief Check a mesh iterator
          * \details  This test performs a series of simple tests on a single mesh element iterator
          * \param[out] ut           Unit test class to report the results
          * \param[in] mesh          mesh for the iterator
          * \param[in] iterator      local iterator over elements
          * \param[in] N_local       number of local elements for the iterator
          * \param[in] N_ghost       number of ghost elements for the iterator
          * \param[in] type          Geometric type
          */
    static void ElementIteratorTest( AMP::UnitTest *ut,
                                     AMP::Mesh::Mesh::shared_ptr mesh,
                                     const AMP::Mesh::MeshIterator &iterator,
                                     const size_t N_local,
                                     const size_t N_ghost,
                                     const AMP::Mesh::GeomType type );

    // Helper function to create a map from the base mesh communicator rank to the main mesh
    // communicator
    static std::map<AMP::Mesh::MeshID, std::vector<int>>
    createRankMap( AMP::Mesh::Mesh::shared_ptr mesh );

    static AMP::Mesh::Mesh::shared_ptr globalMeshForMeshVectorFactory;
};


} // namespace Mesh
} // namespace AMP

// Extra includes
#ifdef USE_AMP_VECTORS
#include "AMP/ampmesh/testHelpers/meshVectorTests.inline.h"
#endif
#ifdef USE_AMP_MATRICES
#include "AMP/ampmesh/testHelpers/meshMatrixTests.inline.h"
#endif


#endif
