#include "DiffusionNonlinearElement.h"
#include "AMP/utils/Utilities.h"
#include <map>
#include <vector>

namespace AMP {
namespace Operator {

void DiffusionNonlinearElement::initTransportModel()
{
    d_fe->reinit( d_elem );

    d_transportModel->preNonlinearInitElementOperation();

    for ( size_t qp = 0; qp < d_qrule->n_points(); qp++ ) {
        d_transportModel->nonlinearInitGaussPointOperation();
    } // end for qp

    d_transportModel->postNonlinearInitElementOperation();
}

void DiffusionNonlinearElement::apply()
{
    const std::vector<Real> &JxW = ( *d_JxW );

    const std::vector<std::vector<Real>> &phi = ( *d_phi );

    const std::vector<std::vector<RealGradient>> &dphi = ( *d_dphi );

    d_fe->reinit( d_elem );

    const size_t num_nodes = d_elem->n_nodes();

    // create transport coefficient storage
    std::vector<double> transportCoeff( d_qrule->n_points() );
    std::vector<std::vector<AMP::shared_ptr<std::vector<double>>>> transportCoeffTensor(
        3, std::vector<AMP::shared_ptr<std::vector<double>>>( 3 ) );
    if ( d_transportModel->isaTensor() ) {
        d_transportTensorModel =
            AMP::dynamic_pointer_cast<DiffusionTransportTensorModel>( d_transportModel );
        for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 3; j++ ) {
                std::vector<double> *vec = new std::vector<double>( d_qrule->n_points() );
                transportCoeffTensor[i][j].reset( vec );
            }
    }

    // compute transport coefficients
    std::map<std::string, AMP::shared_ptr<std::vector<double>>> transport_args;

    // at gauss points
    if ( d_transportAtGauss ) {
        // construct material evalv arguments
        const std::vector<Point> &q_point = d_fe->get_xyz();
        for ( size_t var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
            if ( d_elementInputVectors[var].size() > 0 ) {
                AMP::shared_ptr<std::vector<double>> values(
                    new std::vector<double>( d_qrule->n_points(), 0.0 ) );
                for ( size_t qp = 0; qp < d_qrule->n_points(); qp++ ) {
                    ( *values )[qp] = 0.0;
                    for ( size_t j = 0; j < num_nodes; j++ )
                        ( *values )[qp] += d_elementInputVectors[var][j] * phi[j][qp];
                } // end for qp
                transport_args.insert( std::make_pair( Diffusion::names[var], values ) );
            }
        }

        // evaluate for scalars or tensors
        if ( not d_transportModel->isaTensor() ) {
            d_transportModel->getTransport( transportCoeff, transport_args, q_point );
        } else {
            d_transportTensorModel->getTensorTransport(
                transportCoeffTensor, transport_args, q_point );
        }

        // at nodes
    } else {
        // get element nodes
        std::vector<Point> elem_nodes( num_nodes );
        for ( size_t i = 0; i < num_nodes; i++ ) {
            elem_nodes[i] = d_elem->point( i );
        }

        // set up storage for transport coefficients
        std::vector<std::vector<AMP::shared_ptr<std::vector<double>>>> nodalTransportCoeffTensor(
            3, std::vector<AMP::shared_ptr<std::vector<double>>>( 3 ) );

        // construct material evalv arguments
        for ( size_t var = 0; var < Diffusion::NUMBER_VARIABLES; var++ ) {
            if ( d_elementInputVectors[var].size() > 0 ) {
                transport_args.insert( std::make_pair(
                    Diffusion::names[var],
                    AMP::make_shared<std::vector<double>>( d_elementInputVectors[var] ) ) );
            }
        }

        // evaluate for scalars
        if ( not d_transportModel->isaTensor() ) {
            std::vector<double> nodalTransportCoeff( num_nodes );
            d_transportModel->getTransport( nodalTransportCoeff, transport_args, elem_nodes );

            // interpolate to gauss points
            for ( size_t qp = 0; qp < d_qrule->n_points(); qp++ ) {
                transportCoeff[qp] = 0.0;
                for ( size_t n = 0; n < num_nodes; n++ ) {
                    transportCoeff[qp] += nodalTransportCoeff[n] * phi[n][qp];
                } // end for n
            }     // end for qp

            // evaluate for tensors
        } else {
            for ( int i = 0; i < 3; i++ )
                for ( int j = 0; j < 3; j++ ) {
                    std::vector<double> *vec( new std::vector<double>( num_nodes ) );
                    nodalTransportCoeffTensor[i][j].reset( vec );
                }
            d_transportTensorModel->getTensorTransport(
                nodalTransportCoeffTensor, transport_args, elem_nodes );

            // interpolate to gauss points
            for ( size_t qp = 0; qp < d_qrule->n_points(); qp++ ) {
                for ( int i = 0; i < 3; i++ )
                    for ( int j = 0; j < 3; j++ )
                        ( *transportCoeffTensor[i][j] )[qp] = 0.0;
                for ( size_t n = 0; n < num_nodes; n++ ) {
                    for ( int i = 0; i < 3; i++ )
                        for ( int j = 0; j < 3; j++ ) {
                            ( *transportCoeffTensor[i][j] )[qp] +=
                                ( *nodalTransportCoeffTensor[i][j] )[n] * phi[n][qp];
                        }
                } // end for j
            }     // end for qp
        }
    }

    for ( size_t qp = 0; qp < d_qrule->n_points(); qp++ ) {
        RealGradient grad_phi = 0.0;

        for ( size_t n = 0; n < num_nodes; n++ ) {
            grad_phi += dphi[n][qp] * d_elementInputVectors[d_PrincipalVariable][n];
        }

        d_transportModel->preNonlinearAssemblyGaussPointOperation();

        if ( not d_transportModel->isaTensor() ) {
            for ( size_t n = 0; n < num_nodes; n++ ) {
                ( *d_elementOutputVector )[n] +=
                    ( JxW[qp] * transportCoeff[qp] * dphi[n][qp] * grad_phi );
            } // end for n
        } else {
            for ( size_t n = 0; n < num_nodes; n++ ) {
                for ( int i = 0; i < 3; i++ )
                    for ( int j = 0; j < 3; j++ ) {
                        ( *d_elementOutputVector )[n] +=
                            ( JxW[qp] * ( *transportCoeffTensor[i][j] )[qp] * dphi[n][qp]( i ) *
                              grad_phi( j ) );
                    }
            } // end for n
        }

        d_transportModel->postNonlinearAssemblyGaussPointOperation();
    } // end for qp

    d_transportModel->postNonlinearAssemblyElementOperation();
}
} // namespace Operator
} // namespace AMP
