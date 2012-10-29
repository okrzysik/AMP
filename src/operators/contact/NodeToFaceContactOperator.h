
#ifndef included_AMP_NoteToFaceContactOperator
#define included_AMP_NoteToFaceContactOperator


#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <operators/contact/ContactOperator.h>
#include <matrices/Matrix.h>
#include <vectors/Vector.h>
#include <vectors/Variable.h>

namespace AMP {
  namespace Operator {

    /**
      An abstract base class for representing a linear operator. This class 
      stores the matrix representation of the linear operator. It provides
      an implementation of the apply() function.
      @see Operator
      */
    class NodeToFaceContactOperator : public ContactOperator 
    {

      public :

        /**
          Constructor. This resets the matrix shared pointer.
          @param [in] params 
          */
        NodeToFaceContactOperator (const boost::shared_ptr<ContactOperatorParameters> & params)
          : ContactOperator(params)
        {
          size_t rank = d_GlobalComm.getRank();
          std::string fileName = "debug_operator_" + boost::lexical_cast<std::string>(rank);
          d_fout.open(fileName.c_str(), std::fstream::out);
        }

        /**
          Destructor
          */
        ~NodeToFaceContactOperator() { d_fout.close(); }

        /**
         * This function is useful for re-initializing/updating an operator
         * \param params
         *        parameter object containing parameters to change
         */
        void reset(const boost::shared_ptr<OperatorParameters> & params);

        void addSlaveToMaster(AMP::LinearAlgebra::Vector::shared_ptr u);

        void copyMasterToSlave(AMP::LinearAlgebra::Vector::shared_ptr u); 

        void initialize();

        size_t updateActiveSet();

      protected :

      private :
        void getVectorIndicesFromGlobalIDs(const std::vector<AMP::Mesh::MeshElementID> & globalIDs, 
            std::vector<size_t> & vectorIndices);

        std::vector<int> d_SendCnts;
        std::vector<int> d_SendDisps;
        std::vector<int> d_RecvCnts;
        std::vector<int> d_RecvDisps;
        std::vector<int> d_TransposeSendCnts;
        std::vector<int> d_TransposeSendDisps;
        std::vector<int> d_TransposeRecvCnts;
        std::vector<int> d_TransposeRecvDisps;

        // actually we don't need to store the meshelementids but maybe useful later to check whether the active set has changed
//        std::vector<AMP::Mesh::MeshElementID> d_SlaveVerticesGlobalIDs;
        std::vector<AMP::Mesh::MeshElementID> d_RecvMasterVerticesGlobalIDs;

//        std::vector<size_t> d_SlaveIndices;
        std::vector<size_t> d_RecvMasterIndices;
        std::vector<size_t> d_MasterVerticesMap;
        std::vector<size_t> d_MasterVerticesInverseMap;
        std::vector<size_t> d_MasterVerticesOwnerRanks;
        std::vector<double> d_MasterShapeFunctionsValues;
        std::vector<double> d_SlaveVerticesShift;

        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_InputVariable; /**< Input variable */
        boost::shared_ptr<AMP::LinearAlgebra::Variable> d_OutputVariable; /**< Output variable */

        std::fstream d_fout;
    };

  }
}

#endif


