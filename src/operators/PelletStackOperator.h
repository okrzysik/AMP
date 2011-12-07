
#ifndef included_AMP_PelletStackOperator
#define included_AMP_PelletStackOperator

#include "operators/Operator.h"

namespace AMP {
  namespace Operator {

    class PelletStackOperator : public Operator 
    {
      public :
        PelletStackOperator(const boost::shared_ptr<OperatorParameters> & params);

        ~PelletStackOperator() { }

        bool hasPellet(unsigned int pellId); 

        void setCurrentPellet(unsigned int pellId);

        unsigned int getTotalNumberOfPellets();

        bool useSerial();

        bool onlyZcorrection();

        bool useScaling();

        void applyUnscaling(AMP::LinearAlgebra::Vector::shared_ptr f);

        void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1);

        void setLocalMeshes(std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> inp);

        void setLocalPelletIds(std::vector<unsigned int> inp);

      protected:
        void applySerial(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr  &r);

        void applyOnlyZcorrection(AMP::LinearAlgebra::Vector::shared_ptr &r);

        void applyXYZcorrection(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr  &r);

        unsigned int d_totalNumberOfPellets;
        unsigned int d_currentPellet;
        bool d_useSerial;
        bool d_onlyZcorrection;
        bool d_useScaling;
        double d_scalingFactor;
        short int d_masterId;
        short int d_slaveId;
        std::vector<AMP::Mesh::MeshManager::Adapter::shared_ptr> d_meshes;
        std::vector<unsigned int> d_pelletIds;

    };

  }
}

#endif


