#include "AMP/operators/flow/NavierStokesLSWFLinearElement.h"
#include "AMP/utils/Utilities.h"


namespace AMP::Operator {

void NavierStokesLSWFLinearElement::apply()
{
    const std::vector<libMesh::Real> &JxW = ( *d_JxW );

    const std::vector<std::vector<libMesh::RealGradient>> &dphi = ( *d_dphi );

    const std::vector<std::vector<libMesh::Real>> &phi = ( *d_phi );

    std::vector<std::vector<double>> &elementStiffnessMatrix = ( *d_elementStiffnessMatrix );

    d_fe->reinit( d_elem );

    const unsigned int num_nodes = d_elem->n_nodes();

    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {

        double u = 0;
        double v = 0;
        double w = 0;

        /* double txx = 0;
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
  */
        double dudx = 0;
        double dudy = 0;
        double dudz = 0;
        double dvdx = 0;
        double dvdy = 0;
        double dvdz = 0;
        double dwdx = 0;
        double dwdy = 0;
        double dwdz = 0;

        /* double p = 0;
        double dpdx = 0;
        double dpdy = 0;
        double dpdz = 0;
  */
        for ( unsigned int k = 0; k < num_nodes; k++ ) {

            u += d_elementInputVectors[10 * k + 1] * phi[k][qp];
            v += d_elementInputVectors[10 * k + 2] * phi[k][qp];
            w += d_elementInputVectors[10 * k + 3] * phi[k][qp];
            // p += d_elementInputVectors[ 10*k + 0]*phi[k][qp];


            dudx += ( d_elementInputVectors[( 10 * k ) + 1] * dphi[k][qp]( 0 ) );
            dudy += ( d_elementInputVectors[( 10 * k ) + 1] * dphi[k][qp]( 1 ) );
            dudz += ( d_elementInputVectors[( 10 * k ) + 1] * dphi[k][qp]( 2 ) );

            dvdx += ( d_elementInputVectors[( 10 * k ) + 2] * dphi[k][qp]( 0 ) );
            dvdy += ( d_elementInputVectors[( 10 * k ) + 2] * dphi[k][qp]( 1 ) );
            dvdz += ( d_elementInputVectors[( 10 * k ) + 2] * dphi[k][qp]( 2 ) );

            dwdx += ( d_elementInputVectors[( 10 * k ) + 3] * dphi[k][qp]( 0 ) );
            dwdy += ( d_elementInputVectors[( 10 * k ) + 3] * dphi[k][qp]( 1 ) );
            dwdz += ( d_elementInputVectors[( 10 * k ) + 3] * dphi[k][qp]( 2 ) );

            /*
            dpdx += (d_elementInputVectors[10*k + 0 ]*dphi[k][qp](0));
            dpdy += (d_elementInputVectors[10*k + 0 ]*dphi[k][qp](1));
            dpdz += (d_elementInputVectors[10*k + 0 ]*dphi[k][qp](2));

            txx += d_elementInputVectors[(10*k) + 4]*phi[k][qp];
            tyy += d_elementInputVectors[(10*k) + 5]*phi[k][qp];
            tzz += d_elementInputVectors[(10*k) + 6]*phi[k][qp];

            dtxxdx += (d_elementInputVectors[(10*k) + 4]*dphi[k][qp](0));
            dtyydy += (d_elementInputVectors[(10*k) + 5]*dphi[k][qp](1));
            dtzzdz += (d_elementInputVectors[(10*k) + 6]*dphi[k][qp](2));

            txy += d_elementInputVectors[(10*k) + 7]*phi[k][qp];
            tyz += d_elementInputVectors[(10*k) + 8]*phi[k][qp];
            txz += d_elementInputVectors[(10*k) + 9]*phi[k][qp];

            dtxydx += (d_elementInputVectors[(10*k) + 7]*dphi[k][qp](0));
            dtyzdy += (d_elementInputVectors[(10*k) + 8]*dphi[k][qp](1));
            dtxzdz += (d_elementInputVectors[(10*k) + 9]*dphi[k][qp](2));

            dtxydy += (d_elementInputVectors[(10*k) + 7]*dphi[k][qp](1));
            dtyzdz += (d_elementInputVectors[(10*k) + 8]*dphi[k][qp](2));
            dtxzdx += (d_elementInputVectors[(10*k) + 9]*dphi[k][qp](0));


               u += d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*phi[k][qp];
               v += d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*phi[k][qp];
               w += d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*phi[k][qp];
               p += d_elementInputVectors[NavierStokes::PRESSURE][k]*phi[k][qp];

               dudx += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*dphi[k][qp](0));
               dudy += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*dphi[k][qp](1));
               dudz += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 0]*dphi[k][qp](2));

               dvdx += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*dphi[k][qp](0));
               dvdy += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*dphi[k][qp](1));
               dvdz += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 1]*dphi[k][qp](2));

               dwdx += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*dphi[k][qp](0));
               dwdy += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*dphi[k][qp](1));
               dwdz += (d_elementInputVectors[NavierStokes::VELOCITY][(3*k) + 2]*dphi[k][qp](2));

               dpdx += (d_elementInputVectors[NavierStokes::PRESSURE][k ]*dphi[k][qp](0));
               dpdy += (d_elementInputVectors[NavierStokes::PRESSURE][k ]*dphi[k][qp](1));
               dpdz += (d_elementInputVectors[NavierStokes::PRESSURE][k ]*dphi[k][qp](2));

               txx += d_elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 0]*phi[k][qp];
               tyy += d_elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 1]*phi[k][qp];
               tzz += d_elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) + 2]*phi[k][qp];

               dtxxdx += (d_elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) +
            0]*dphi[k][qp](0));
               dtyydy += (d_elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) +
            1]*dphi[k][qp](1));
               dtzzdz += (d_elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*k) +
            2]*dphi[k][qp](2));

               txy += d_elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 0]*phi[k][qp];
               tyz += d_elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 1]*phi[k][qp];
               txz += d_elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) + 2]*phi[k][qp];

               dtxydx += (d_elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
            0]*dphi[k][qp](0));
               dtyzdy += (d_elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
            1]*dphi[k][qp](1));
               dtxzdz += (d_elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
            2]*dphi[k][qp](2));

               dtxydy += (d_elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
            0]*dphi[k][qp](1));
               dtyzdz += (d_elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
            1]*dphi[k][qp](2));
               dtxzdx += (d_elementInputVectors[NavierStokes::SHEARSTRESS][(3*k) +
            2]*dphi[k][qp](0));
               */

        } // end for k

        d_density = d_transportModel->getDensity();
        d_fmu     = d_transportModel->getViscosity();
        d_Re      = d_transportModel->getReynoldsNumber();

        std::vector<std::vector<double>> derror_i( 10, std::vector<double>( 10, 0 ) );
        std::vector<std::vector<double>> derror_j( 10, std::vector<double>( 10, 0 ) );

        for ( unsigned int i = 0; i < num_nodes; i++ ) {

            derror_i[0][1] = dphi[i][qp]( 0 );
            derror_i[0][2] = dphi[i][qp]( 1 );
            derror_i[0][3] = dphi[i][qp]( 2 );

            derror_i[1][0] = dphi[i][qp]( 0 );
            derror_i[1][1] = ( phi[i][qp] * dudx + u * dphi[i][qp]( 0 ) + v * dphi[i][qp]( 1 ) +
                               w * dphi[i][qp]( 2 ) );
            derror_i[1][2] = ( phi[i][qp] * dudy );
            derror_i[1][3] = ( phi[i][qp] * dudz );
            derror_i[1][4] = -1 * dphi[i][qp]( 0 );
            derror_i[1][7] = -1 * dphi[i][qp]( 1 );
            derror_i[1][9] = -1 * dphi[i][qp]( 2 );

            derror_i[2][0] = dphi[i][qp]( 1 );
            derror_i[2][1] = ( phi[i][qp] * dvdx );
            derror_i[2][2] = ( phi[i][qp] * dvdy + v * dphi[i][qp]( 1 ) + u * dphi[i][qp]( 0 ) +
                               w * dphi[i][qp]( 2 ) );
            derror_i[2][3] = ( phi[i][qp] * dvdz );
            derror_i[2][5] = -1 * dphi[i][qp]( 1 );
            derror_i[2][7] = -1 * dphi[i][qp]( 0 );
            derror_i[2][8] = -1 * dphi[i][qp]( 2 );

            derror_i[3][0] = dphi[i][qp]( 2 );
            derror_i[3][1] = ( phi[i][qp] * dwdx );
            derror_i[3][2] = ( phi[i][qp] * dwdy );
            derror_i[3][3] = ( phi[i][qp] * dwdz + w * dphi[i][qp]( 2 ) + u * dphi[i][qp]( 0 ) +
                               v * dphi[i][qp]( 1 ) );
            derror_i[3][6] = -1 * dphi[i][qp]( 2 );
            derror_i[3][8] = -1 * dphi[i][qp]( 1 );
            derror_i[3][9] = -1 * dphi[i][qp]( 0 );

            derror_i[4][1] = -1 * ( 2.0 * ( d_fmu / d_Re ) * dphi[i][qp]( 0 ) );
            derror_i[4][4] = phi[i][qp];

            derror_i[5][2] = -1 * ( 2.0 * ( d_fmu / d_Re ) * dphi[i][qp]( 1 ) );
            derror_i[5][5] = phi[i][qp];

            derror_i[6][3] = -1 * ( 2.0 * ( d_fmu / d_Re ) * dphi[i][qp]( 2 ) );
            derror_i[6][6] = phi[i][qp];

            derror_i[7][1] = -1 * ( ( d_fmu / d_Re ) * dphi[i][qp]( 1 ) );
            derror_i[7][2] = -1 * ( ( d_fmu / d_Re ) * dphi[i][qp]( 0 ) );
            derror_i[7][7] = phi[i][qp];

            derror_i[8][2] = -1 * ( ( d_fmu / d_Re ) * dphi[i][qp]( 2 ) );
            derror_i[8][3] = -1 * ( ( d_fmu / d_Re ) * dphi[i][qp]( 1 ) );
            derror_i[8][8] = phi[i][qp];

            derror_i[9][1] = -1 * ( ( d_fmu / d_Re ) * dphi[i][qp]( 2 ) );
            derror_i[9][3] = -1 * ( ( d_fmu / d_Re ) * dphi[i][qp]( 0 ) );
            derror_i[9][9] = phi[i][qp];

            for ( unsigned int j = 0; j < num_nodes; j++ ) {

                derror_j[0][1] = dphi[j][qp]( 0 );
                derror_j[0][2] = dphi[j][qp]( 1 );
                derror_j[0][3] = dphi[j][qp]( 2 );

                derror_j[1][0] = dphi[j][qp]( 0 );
                derror_j[1][1] = ( phi[j][qp] * dudx + u * dphi[j][qp]( 0 ) + v * dphi[j][qp]( 1 ) +
                                   w * dphi[j][qp]( 2 ) );
                derror_j[1][2] = ( phi[j][qp] * dudy );
                derror_j[1][3] = ( phi[j][qp] * dudz );
                derror_j[1][4] = -1 * dphi[j][qp]( 0 );
                derror_j[1][7] = -1 * dphi[j][qp]( 1 );
                derror_j[1][9] = -1 * dphi[j][qp]( 2 );

                derror_j[2][0] = dphi[j][qp]( 1 );
                derror_j[2][1] = ( phi[j][qp] * dvdx );
                derror_j[2][2] = ( phi[j][qp] * dvdy + v * dphi[j][qp]( 1 ) + u * dphi[j][qp]( 0 ) +
                                   w * dphi[j][qp]( 2 ) );
                derror_j[2][3] = ( phi[j][qp] * dvdz );
                derror_j[2][5] = -1 * dphi[j][qp]( 1 );
                derror_j[2][7] = -1 * dphi[j][qp]( 0 );
                derror_j[2][8] = -1 * dphi[j][qp]( 2 );

                derror_j[3][0] = dphi[j][qp]( 2 );
                derror_j[3][1] = ( phi[j][qp] * dwdx );
                derror_j[3][2] = ( phi[j][qp] * dwdy );
                derror_j[3][3] = ( phi[j][qp] * dwdz + w * dphi[j][qp]( 2 ) + u * dphi[j][qp]( 0 ) +
                                   v * dphi[j][qp]( 1 ) );
                derror_j[3][6] = -1 * dphi[j][qp]( 2 );
                derror_j[3][8] = -1 * dphi[j][qp]( 1 );
                derror_j[3][9] = -1 * dphi[j][qp]( 0 );

                derror_j[4][1] = -1 * ( 2.0 * ( d_fmu / d_Re ) * dphi[j][qp]( 0 ) );
                derror_j[4][4] = phi[j][qp];

                derror_j[5][2] = -1 * ( 2.0 * ( d_fmu / d_Re ) * dphi[j][qp]( 1 ) );
                derror_j[5][5] = phi[j][qp];

                derror_j[6][3] = -1 * ( 2.0 * ( d_fmu / d_Re ) * dphi[j][qp]( 2 ) );
                derror_j[6][6] = phi[j][qp];

                derror_j[7][1] = -1 * ( ( d_fmu / d_Re ) * dphi[j][qp]( 1 ) );
                derror_j[7][2] = -1 * ( ( d_fmu / d_Re ) * dphi[j][qp]( 0 ) );
                derror_j[7][7] = phi[j][qp];

                derror_j[8][2] = -1 * ( ( d_fmu / d_Re ) * dphi[j][qp]( 2 ) );
                derror_j[8][3] = -1 * ( ( d_fmu / d_Re ) * dphi[j][qp]( 1 ) );
                derror_j[8][8] = phi[j][qp];

                derror_j[9][1] = -1 * ( ( d_fmu / d_Re ) * dphi[j][qp]( 2 ) );
                derror_j[9][3] = -1 * ( ( d_fmu / d_Re ) * dphi[j][qp]( 0 ) );
                derror_j[9][9] = phi[j][qp];

                for ( unsigned int ii = 0; ii < 10; ii++ ) {

                    for ( unsigned int jj = 0; jj < 10; jj++ ) {

                        for ( unsigned int eq = 0; eq < 10; eq++ )
                            elementStiffnessMatrix[10 * i + ii][10 * j + jj] +=
                                JxW[qp] * ( derror_i[eq][ii] * derror_j[eq][jj] );

                    } // jj
                }     // ii

            } // end for j
        }     // end for i

    } // end for qp
}
} // namespace AMP::Operator
