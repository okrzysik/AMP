
#ifndef included_AMP_PelletStackOperator
#define included_AMP_PelletStackOperator

#include "AMP/operators/libmesh/PelletStackOperatorParameters.h"

namespace AMP {
namespace Operator {

class PelletStackOperator : public Operator
{
public:
    explicit PelletStackOperator( std::shared_ptr<const PelletStackOperatorParameters> params );

    virtual ~PelletStackOperator() {}

    std::string type() const override { return "PelletStackOperator"; }

    int getLocalIndexForPellet( unsigned int pellId );

    unsigned int getTotalNumberOfPellets();

    std::vector<AMP::Mesh::Mesh::shared_ptr> getLocalMeshes();

    std::vector<unsigned int> getLocalPelletIds();

    bool useSerial();

    bool onlyZcorrection();

    bool useScaling();

    void reset( std::shared_ptr<const OperatorParameters> params ) override;

    void applyUnscaling( AMP::LinearAlgebra::Vector::shared_ptr f );

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

protected:
    void applySerial( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                      AMP::LinearAlgebra::Vector::shared_ptr &f );

    void applyOnlyZcorrection( AMP::LinearAlgebra::Vector::shared_ptr &u );

    void applyXYZcorrection( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr &f );

    void computeZscan( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                       std::vector<double> &finalMaxZdispsList );

    unsigned int d_totalNumberOfPellets;
    unsigned int d_currentPellet;
    bool d_useSerial;
    bool d_onlyZcorrection;
    bool d_useScaling;
    double d_scalingFactor;
    short int d_masterId;
    short int d_slaveId;
    std::vector<AMP::Mesh::Mesh::shared_ptr> d_meshes;
    std::vector<unsigned int> d_pelletIds;
    AMP::LinearAlgebra::Variable::shared_ptr d_var;
    AMP::LinearAlgebra::Vector::shared_ptr d_frozenVectorForMaps;
    bool d_frozenVectorSet;
    AMP_MPI d_pelletStackComm;
    std::shared_ptr<AMP::Operator::AsyncMapColumnOperator> d_n2nMaps;
};
} // namespace Operator
} // namespace AMP

#endif
