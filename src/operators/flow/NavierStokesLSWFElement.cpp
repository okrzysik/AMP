
#include "NavierStokesLSWFElement.h"
#include "AMP/utils/Utilities.h"

namespace AMP::Operator {

void NavierStokesLSWFElement::initTransportModel() { d_fe->reinit( d_elem ); }

void NavierStokesLSWFElement::apply()
{
    const std::vector<libMesh::Real> &JxW = ( *d_JxW );

    const std::vector<std::vector<libMesh::RealGradient>> &dphi = ( *d_dphi );

    const std::vector<std::vector<libMesh::Real>> &phi = ( *d_phi );

    std::vector<double> &elementInputVectors = d_elementInputVectors;

    std::vector<double> &elementOutputVector = *d_elementOutputVector;

    d_fe->reinit( d_elem );

    d_fmu = d_transportModel->getViscosity();

    double d_Re = d_transportModel->getReynoldsNumber();
    NULL_USE( d_Re );

    const unsigned int num_nodes = d_elem->n_nodes();

    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {

        double u = 0;
        double v = 0;
        double w = 0;

        double txx = 0;
        double tyy = 0;
        double tzz = 0;
        double txy = 0;
        double tyz = 0;
        double txz = 0;

        double dtxxdx = 0;
        double dtyydy = 0;
        double dtzzdz = 0;

        double dtxydx = 0;
        double dtxydy = 0;

        double dtyzdy = 0;
        double dtyzdz = 0;

        double dtxzdx = 0;
        double dtxzdz = 0;

        double dudx = 0;
        double dudy = 0;
        double dudz = 0;
        double dvdx = 0;
        double dvdy = 0;
        double dvdz = 0;
        double dwdx = 0;
        double dwdy = 0;
        double dwdz = 0;

        // double p = 0;
        double dpdx = 0;
        double dpdy = 0;
        double dpdz = 0;

        for ( unsigned int k = 0; k < num_nodes; k++ ) {
            u += elementInputVectors[10 * k + 1] * phi[k][qp];
            v += elementInputVectors[10 * k + 2] * phi[k][qp];
            w += elementInputVectors[10 * k + 3] * phi[k][qp];
            // p += elementInputVectors[ 10*k + 0]*phi[k][qp]; // code that uses it is commented
            // out.


            dudx += ( elementInputVectors[( 10 * k ) + 1] * dphi[k][qp]( 0 ) );
            dudy += ( elementInputVectors[( 10 * k ) + 1] * dphi[k][qp]( 1 ) );
            dudz += ( elementInputVectors[( 10 * k ) + 1] * dphi[k][qp]( 2 ) );

            dvdx += ( elementInputVectors[( 10 * k ) + 2] * dphi[k][qp]( 0 ) );
            dvdy += ( elementInputVectors[( 10 * k ) + 2] * dphi[k][qp]( 1 ) );
            dvdz += ( elementInputVectors[( 10 * k ) + 2] * dphi[k][qp]( 2 ) );

            dwdx += ( elementInputVectors[( 10 * k ) + 3] * dphi[k][qp]( 0 ) );
            dwdy += ( elementInputVectors[( 10 * k ) + 3] * dphi[k][qp]( 1 ) );
            dwdz += ( elementInputVectors[( 10 * k ) + 3] * dphi[k][qp]( 2 ) );

            dpdx += ( elementInputVectors[10 * k + 0] * dphi[k][qp]( 0 ) );
            dpdy += ( elementInputVectors[10 * k + 0] * dphi[k][qp]( 1 ) );
            dpdz += ( elementInputVectors[10 * k + 0] * dphi[k][qp]( 2 ) );

            txx += elementInputVectors[( 10 * k ) + 4] * phi[k][qp];
            tyy += elementInputVectors[( 10 * k ) + 5] * phi[k][qp];
            tzz += elementInputVectors[( 10 * k ) + 6] * phi[k][qp];

            dtxxdx += ( elementInputVectors[( 10 * k ) + 4] * dphi[k][qp]( 0 ) );
            dtyydy += ( elementInputVectors[( 10 * k ) + 5] * dphi[k][qp]( 1 ) );
            dtzzdz += ( elementInputVectors[( 10 * k ) + 6] * dphi[k][qp]( 2 ) );

            txy += elementInputVectors[( 10 * k ) + 7] * phi[k][qp];
            tyz += elementInputVectors[( 10 * k ) + 8] * phi[k][qp];
            txz += elementInputVectors[( 10 * k ) + 9] * phi[k][qp];

            dtxydx += ( elementInputVectors[( 10 * k ) + 7] * dphi[k][qp]( 0 ) );
            dtyzdy += ( elementInputVectors[( 10 * k ) + 8] * dphi[k][qp]( 1 ) );
            dtxzdz += ( elementInputVectors[( 10 * k ) + 9] * dphi[k][qp]( 2 ) );

            dtxydy += ( elementInputVectors[( 10 * k ) + 7] * dphi[k][qp]( 1 ) );
            dtyzdz += ( elementInputVectors[( 10 * k ) + 8] * dphi[k][qp]( 2 ) );
            dtxzdx += ( elementInputVectors[( 10 * k ) + 9] * dphi[k][qp]( 0 ) );


            /*
                    u += elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*phi[k][qp];
                    v += elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*phi[k][qp];
                    w += elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*phi[k][qp];
                    p += elementInputVectors[NavierStokes::PRESSURE][k]*phi[k][qp];


                    dudx += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*dphi[k][qp](0));
                    dudy += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*dphi[k][qp](1));
                    dudz += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*dphi[k][qp](2));

                    dvdx += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*dphi[k][qp](0));
                    dvdy += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*dphi[k][qp](1));
                    dvdz += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*dphi[k][qp](2));

                    dwdx += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*dphi[k][qp](0));
                    dwdy += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*dphi[k][qp](1));
                    dwdz += (elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*dphi[k][qp](2));

                    dpdx += (elementInputVectors[NavierStokes::PRESSURE][k ]*dphi[k][qp](0));
                    dpdy += (elementInputVectors[NavierStokes::PRESSURE][k ]*dphi[k][qp](1));
                    dpdz += (elementInputVectors[NavierStokes::PRESSURE][k ]*dphi[k][qp](2));

                    txx += elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 0]*phi[k][qp];
                    tyy += elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 1]*phi[k][qp];
                    tzz += elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 2]*phi[k][qp];

                    dtxxdx += (elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) +
               0]*dphi[k][qp](0));
                    dtyydy += (elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) +
               1]*dphi[k][qp](1));
                    dtzzdz += (elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) +
               2]*dphi[k][qp](2));

                    txy += elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 0]*phi[k][qp];
                    tyz += elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 1]*phi[k][qp];
                    txz += elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 2]*phi[k][qp];

                    dtxydx += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
               0]*dphi[k][qp](0));
                    dtyzdy += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
               1]*dphi[k][qp](1));
                    dtxzdz += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
               2]*dphi[k][qp](2));

                    dtxydy += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
               0]*dphi[k][qp](1));
                    dtyzdz += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
               1]*dphi[k][qp](2));
                    dtxzdx += (elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
               2]*dphi[k][qp](0));
            */
        } // end for k

        d_density = d_transportModel->getDensity();
        d_fmu     = d_transportModel->getViscosity();
        d_Re      = d_transportModel->getReynoldsNumber();

        std::vector<double> error( 10 );
        std::vector<std::vector<double>> derror( 10, std::vector<double>( 10, 0 ) );

        error[0] = ( dudx + dvdy + dwdz );
        error[1] = ( u * dudx + v * dudy + w * dudz ) + dpdx - ( dtxxdx + dtxydy + dtxzdz );
        error[2] = ( u * dvdx + v * dvdy + w * dvdz ) + dpdy - ( dtxydx + dtyydy + dtyzdz );
        error[3] = ( u * dwdx + v * dwdy + w * dwdz ) + dpdz - ( dtxzdx + dtyzdy + dtzzdz );
        error[4] = txx - 2.0 * ( d_fmu / d_Re ) * dudx;
        error[5] = tyy - 2.0 * ( d_fmu / d_Re ) * dvdy;
        error[6] = tzz - 2.0 * ( d_fmu / d_Re ) * dwdz;
        error[7] = txy - ( d_fmu / d_Re ) * ( dudy + dvdx );
        error[8] = tyz - ( d_fmu / d_Re ) * ( dvdz + dwdy );
        error[9] = txz - ( d_fmu / d_Re ) * ( dwdx + dudz );

        for ( unsigned int k = 0; k < num_nodes; k++ ) {
            derror[0][1] = dphi[k][qp]( 0 );
            derror[0][2] = dphi[k][qp]( 1 );
            derror[0][3] = dphi[k][qp]( 2 );

            derror[1][0] = dphi[k][qp]( 0 );
            derror[1][1] = ( phi[k][qp] * dudx + u * dphi[k][qp]( 0 ) + v * dphi[k][qp]( 1 ) +
                             w * dphi[k][qp]( 2 ) );
            derror[1][2] = ( phi[k][qp] * dudy );
            derror[1][3] = ( phi[k][qp] * dudz );
            derror[1][4] = -1 * dphi[k][qp]( 0 );
            derror[1][7] = -1 * dphi[k][qp]( 1 );
            derror[1][9] = -1 * dphi[k][qp]( 2 );

            derror[2][0] = dphi[k][qp]( 1 );
            derror[2][1] = ( phi[k][qp] * dvdx );
            derror[2][2] = ( phi[k][qp] * dvdy + v * dphi[k][qp]( 1 ) + u * dphi[k][qp]( 0 ) +
                             w * dphi[k][qp]( 2 ) );
            derror[2][3] = ( phi[k][qp] * dvdz );
            derror[2][5] = -1 * dphi[k][qp]( 1 );
            derror[2][7] = -1 * dphi[k][qp]( 0 );
            derror[2][8] = -1 * dphi[k][qp]( 2 );

            derror[3][0] = dphi[k][qp]( 2 );
            derror[3][1] = ( phi[k][qp] * dwdx );
            derror[3][2] = ( phi[k][qp] * dwdy );
            derror[3][3] = ( phi[k][qp] * dwdz + w * dphi[k][qp]( 2 ) + u * dphi[k][qp]( 0 ) +
                             v * dphi[k][qp]( 1 ) );
            derror[3][6] = -1 * dphi[k][qp]( 2 );
            derror[3][8] = -1 * dphi[k][qp]( 1 );
            derror[3][9] = -1 * dphi[k][qp]( 0 );

            derror[4][1] = -1 * ( 2.0 * ( d_fmu / d_Re ) * dphi[k][qp]( 0 ) );
            derror[4][4] = phi[k][qp];

            derror[5][2] = -1 * ( 2.0 * ( d_fmu / d_Re ) * dphi[k][qp]( 1 ) );
            derror[5][5] = phi[k][qp];

            derror[6][3] = -1 * ( 2.0 * ( d_fmu / d_Re ) * dphi[k][qp]( 2 ) );
            derror[6][6] = phi[k][qp];

            derror[7][1] = -1 * ( ( d_fmu / d_Re ) * dphi[k][qp]( 1 ) );
            derror[7][2] = -1 * ( ( d_fmu / d_Re ) * dphi[k][qp]( 0 ) );
            derror[7][7] = phi[k][qp];

            derror[8][2] = -1 * ( ( d_fmu / d_Re ) * dphi[k][qp]( 2 ) );
            derror[8][3] = -1 * ( ( d_fmu / d_Re ) * dphi[k][qp]( 1 ) );
            derror[8][8] = phi[k][qp];

            derror[9][1] = -1 * ( ( d_fmu / d_Re ) * dphi[k][qp]( 2 ) );
            derror[9][3] = -1 * ( ( d_fmu / d_Re ) * dphi[k][qp]( 0 ) );
            derror[9][9] = phi[k][qp];

            for ( unsigned int ii = 0; ii < 10; ii++ ) {
                elementOutputVector[10 * k] += JxW[qp] * ( error[0] * derror[0][ii] );
                elementOutputVector[10 * k + 1] += JxW[qp] * ( error[1] * derror[1][ii] );
                elementOutputVector[10 * k + 2] += JxW[qp] * ( error[2] * derror[2][ii] );
                elementOutputVector[10 * k + 3] += JxW[qp] * ( error[3] * derror[3][ii] );
                elementOutputVector[10 * k + 4] += JxW[qp] * ( error[4] * derror[4][ii] );
                elementOutputVector[10 * k + 5] += JxW[qp] * ( error[5] * derror[5][ii] );
                elementOutputVector[10 * k + 6] += JxW[qp] * ( error[6] * derror[6][ii] );
                elementOutputVector[10 * k + 7] += JxW[qp] * ( error[7] * derror[7][ii] );
                elementOutputVector[10 * k + 8] += JxW[qp] * ( error[8] * derror[8][ii] );
                elementOutputVector[10 * k + 9] += JxW[qp] * ( error[9] * derror[9][ii] );
            }
        }
    } // end for qp
}
} // namespace AMP::Operator
