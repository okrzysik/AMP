
#include "NonlinearFEOperator.h"
#include "utils/Utilities.h"
#include "utils/ProfilerApp.h"
#include "cell_hex8.h"
#include "node.h"

namespace AMP {
namespace Operator {

void NonlinearFEOperator :: apply(const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &f, 
    const boost::shared_ptr<AMP::LinearAlgebra::Vector>  &u, boost::shared_ptr<AMP::LinearAlgebra::Vector>  &r,
    const double a,  const double b)
{
    PROFILE_START("apply");
    AMP_INSIST( (r != NULL), "NULL Residual/Output Vector" );

    AMP::LinearAlgebra::Vector::shared_ptr rInternal = this->subsetOutputVector(r);

    AMP_INSIST( (rInternal != NULL), "NULL Residual/Output Vector" );

    AMP::Mesh::MeshIterator el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
    AMP::Mesh::MeshIterator end_el = el.end();

    d_currElemIdx = static_cast<unsigned int>(-1);
    this->preAssembly(u, rInternal);

    PROFILE_START("loop over elements");
    for(d_currElemIdx = 0; el != end_el; ++el, ++d_currElemIdx) {
        this->preElementOperation(*el);
        d_elemOp->apply();
        this->postElementOperation();
    }//end for el
    PROFILE_STOP("loop over elements");

    d_currElemIdx = static_cast<unsigned int>(-1);
    this->postAssembly();

    if(f == NULL) {
        rInternal->scale(a);
    } else {
        AMP::LinearAlgebra::Vector::shared_ptr fInternal = this->subsetOutputVector(f);
        if(fInternal == NULL) {
            rInternal->scale(a);
        } else {
            rInternal->axpby(b, a, fInternal);
        }
    }

    if(d_iDebugPrintInfoLevel>2)
        AMP::pout << "L2 norm of result of NonlinearFEOperator::apply is: " << rInternal->L2Norm() << std::endl;
    if(d_iDebugPrintInfoLevel>5)
        std::cout << rInternal << std::endl;

    PROFILE_STOP("apply");
}


void NonlinearFEOperator :: createLibMeshElementList() {
    AMP::Mesh::MeshIterator el = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
    AMP::Mesh::MeshIterator end_el = el.end();
    d_currElemPtrs.resize(d_Mesh->numLocalElements(AMP::Mesh::Volume));
    for(size_t i = 0; el != end_el; ++el, ++i) {
        std::vector<AMP::Mesh::MeshElement> currNodes = el->getElements(AMP::Mesh::Vertex);
        d_currElemPtrs[i] = new ::Hex8;
        for(size_t j = 0; j < currNodes.size(); ++j) {
            std::vector<double> pt = currNodes[j].coord();
            d_currElemPtrs[i]->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
        }//end for j
    }//end for i
}


void NonlinearFEOperator :: destroyLibMeshElementList() {
    for(size_t i = 0; i < d_currElemPtrs.size(); ++i) {
        for(size_t j = 0; j < d_currElemPtrs[i]->n_nodes(); ++j) {
            delete (d_currElemPtrs[i]->get_node(j));
            d_currElemPtrs[i]->set_node(j) = NULL;
        }//end for j
        delete (d_currElemPtrs[i]);
        d_currElemPtrs[i] = NULL;
    }//end for i
    d_currElemPtrs.clear();
}


}
}//end namespace


