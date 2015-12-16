
#include "MechanicsLinearElement.h"
#include "utils/Utilities.h"

namespace AMP {
namespace Operator {

void MechanicsLinearElement::computeStressAndStrain( const std::vector<double> &dispVec,
                                                     std::vector<double> &stressVec,
                                                     std::vector<double> &strainVec )
{
    const std::vector<std::vector<RealGradient>> &dphi = ( *d_dphi );

    d_fe->reinit( d_elem );

    d_materialModel->preLinearElementOperation();

    const unsigned int num_nodes = d_elem->n_nodes();

    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
        d_materialModel->preLinearGaussPointOperation();

        double *constitutiveMatrix;

        d_materialModel->getConstitutiveMatrix( constitutiveMatrix );

        /* Compute Strain From Given Displacement */

        double dudx = 0;
        double dudy = 0;
        double dudz = 0;
        double dvdx = 0;
        double dvdy = 0;
        double dvdz = 0;
        double dwdx = 0;
        double dwdy = 0;
        double dwdz = 0;

        for ( unsigned int k = 0; k < num_nodes; k++ ) {
            dudx += ( dispVec[( 3 * k ) + 0] * dphi[k][qp]( 0 ) );
            dudy += ( dispVec[( 3 * k ) + 0] * dphi[k][qp]( 1 ) );
            dudz += ( dispVec[( 3 * k ) + 0] * dphi[k][qp]( 2 ) );

            dvdx += ( dispVec[( 3 * k ) + 1] * dphi[k][qp]( 0 ) );
            dvdy += ( dispVec[( 3 * k ) + 1] * dphi[k][qp]( 1 ) );
            dvdz += ( dispVec[( 3 * k ) + 1] * dphi[k][qp]( 2 ) );

            dwdx += ( dispVec[( 3 * k ) + 2] * dphi[k][qp]( 0 ) );
            dwdy += ( dispVec[( 3 * k ) + 2] * dphi[k][qp]( 1 ) );
            dwdz += ( dispVec[( 3 * k ) + 2] * dphi[k][qp]( 2 ) );
        } // end for k

        double uStrain[6];
        uStrain[0] = dudx;
        uStrain[1] = dvdy;
        uStrain[2] = dwdz;
        // uStrain[3] = 0.5*(dvdz + dwdy);
        // uStrain[4] = 0.5*(dudz + dwdx);
        // uStrain[5] = 0.5*(dudy + dvdx);
        uStrain[3] = ( dvdz + dwdy );
        uStrain[4] = ( dudz + dwdx );
        uStrain[5] = ( dudy + dvdx );

        double uStress[6];
        for ( int r = 0; r < 6; r++ ) {
            uStress[r] = 0;
            for ( int c = 0; c < 6; c++ ) {
                uStress[r] += ( constitutiveMatrix[( 6 * r ) + c] * uStrain[c] );
            } // end for c
        }     // end for r

        for ( int i = 0; i < 6; i++ ) {
            stressVec[( 6 * qp ) + i] = uStress[i];
            strainVec[( 6 * qp ) + i] = uStrain[i];
        } // end for i

        d_materialModel->postLinearGaussPointOperation();
    } // end for qp

    d_materialModel->postLinearElementOperation();
}

void MechanicsLinearElement::printStressAndStrain( FILE *fp, const std::vector<double> &dispVec )
{
    const std::vector<std::vector<RealGradient>> &dphi = ( *d_dphi );

    const std::vector<Point> &xyz = ( *d_xyz );

    d_fe->reinit( d_elem );

    d_materialModel->preLinearElementOperation();

    const unsigned int num_nodes = d_elem->n_nodes();

    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
        d_materialModel->preLinearGaussPointOperation();

        double *constitutiveMatrix;

        d_materialModel->getConstitutiveMatrix( constitutiveMatrix );

        /* Compute Strain From Given Displacement */

        double dudx = 0;
        double dudy = 0;
        double dudz = 0;
        double dvdx = 0;
        double dvdy = 0;
        double dvdz = 0;
        double dwdx = 0;
        double dwdy = 0;
        double dwdz = 0;

        for ( unsigned int k = 0; k < num_nodes; k++ ) {
            dudx += ( dispVec[( 3 * k ) + 0] * dphi[k][qp]( 0 ) );
            dudy += ( dispVec[( 3 * k ) + 0] * dphi[k][qp]( 1 ) );
            dudz += ( dispVec[( 3 * k ) + 0] * dphi[k][qp]( 2 ) );

            dvdx += ( dispVec[( 3 * k ) + 1] * dphi[k][qp]( 0 ) );
            dvdy += ( dispVec[( 3 * k ) + 1] * dphi[k][qp]( 1 ) );
            dvdz += ( dispVec[( 3 * k ) + 1] * dphi[k][qp]( 2 ) );

            dwdx += ( dispVec[( 3 * k ) + 2] * dphi[k][qp]( 0 ) );
            dwdy += ( dispVec[( 3 * k ) + 2] * dphi[k][qp]( 1 ) );
            dwdz += ( dispVec[( 3 * k ) + 2] * dphi[k][qp]( 2 ) );
        } // end for k

        double uStrain[6];
        uStrain[0] = dudx;
        uStrain[1] = dvdy;
        uStrain[2] = dwdz;
        // uStrain[3] = 0.5*(dvdz + dwdy);
        // uStrain[4] = 0.5*(dudz + dwdx);
        // uStrain[5] = 0.5*(dudy + dvdx);
        uStrain[3] = ( dvdz + dwdy );
        uStrain[4] = ( dudz + dwdx );
        uStrain[5] = ( dudy + dvdx );

        double uStress[6];
        for ( int r = 0; r < 6; r++ ) {
            uStress[r] = 0;
            for ( int c = 0; c < 6; c++ ) {
                uStress[r] += ( constitutiveMatrix[( 6 * r ) + c] * uStrain[c] );
            } // end for c
        }     // end for r

        fprintf( fp, "%.12f %.12f %.12f \n", xyz[qp]( 0 ), xyz[qp]( 1 ), xyz[qp]( 2 ) );
        fprintf( fp,
                 "%.12f %.12f %.12f %.12f %.12f %.12f \n",
                 uStress[0],
                 uStress[1],
                 uStress[2],
                 uStress[3],
                 uStress[4],
                 uStress[5] );
        fprintf( fp,
                 "%.12f %.12f %.12f %.12f %.12f %.12f \n\n",
                 uStrain[0],
                 uStrain[1],
                 uStrain[2],
                 uStrain[3],
                 uStrain[4],
                 uStrain[5] );

        d_materialModel->postLinearGaussPointOperation();
    } // end for qp

    d_materialModel->postLinearElementOperation();
}

void MechanicsLinearElement::apply_Reduced()
{
    const std::vector<Real> &JxW = ( *d_JxW );

    const std::vector<std::vector<RealGradient>> &dphi = ( *d_dphi );

    std::vector<std::vector<double>> &elementStiffnessMatrix = ( *d_elementStiffnessMatrix );

    d_fe->reinit( d_elem );

    d_materialModel->preLinearElementOperation();

    const unsigned int num_nodes = d_elem->n_nodes();

    std::vector<double> avgTraceTerm( num_nodes * 3 );
    for ( unsigned int i = 0; i < ( 3 * num_nodes ); i++ ) {
        avgTraceTerm[i] = 0.0;
    } // end for i

    double sumJxW = 0.0;
    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
        sumJxW += JxW[qp];

        for ( unsigned int j = 0; j < num_nodes; j++ ) {
            for ( int d = 0; d < 3; d++ ) {
                avgTraceTerm[( 3 * j ) + d] += ( JxW[qp] * dphi[j][qp]( d ) );
            } // end for d
        }     // end for j
    }         // end for qp

    for ( unsigned int i = 0; i < ( 3 * num_nodes ); i++ ) {
        avgTraceTerm[i] = avgTraceTerm[i] / sumJxW;
    } // end for i

    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
        d_materialModel->preLinearGaussPointOperation();

        double *constitutiveMatrix;
        double Bl_np1[6][24], materialMatrix[6][6];
        double materialStiffness[24][24], materialStiffnessTemp[6][24];

        for ( int i = 0; i < 6; i++ ) {
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                Bl_np1[i][j] = 0.0;
            }
        }

        for ( unsigned int i = 0; i < num_nodes; i++ ) {
            double one3 = 1.0 / 3.0;
            Bl_np1[0][( 3 * i ) + 0] =
                dphi[i][qp]( 0 ) + ( one3 * ( avgTraceTerm[( 3 * i ) + 0] - dphi[i][qp]( 0 ) ) );
            Bl_np1[1][( 3 * i ) + 0] = one3 * ( avgTraceTerm[( 3 * i ) + 0] - dphi[i][qp]( 0 ) );
            Bl_np1[2][( 3 * i ) + 0] = one3 * ( avgTraceTerm[( 3 * i ) + 0] - dphi[i][qp]( 0 ) );
            Bl_np1[0][( 3 * i ) + 1] = one3 * ( avgTraceTerm[( 3 * i ) + 1] - dphi[i][qp]( 1 ) );
            Bl_np1[1][( 3 * i ) + 1] =
                dphi[i][qp]( 1 ) + ( one3 * ( avgTraceTerm[( 3 * i ) + 1] - dphi[i][qp]( 1 ) ) );
            Bl_np1[2][( 3 * i ) + 1] = one3 * ( avgTraceTerm[( 3 * i ) + 1] - dphi[i][qp]( 1 ) );
            Bl_np1[0][( 3 * i ) + 2] = one3 * ( avgTraceTerm[( 3 * i ) + 2] - dphi[i][qp]( 2 ) );
            Bl_np1[1][( 3 * i ) + 2] = one3 * ( avgTraceTerm[( 3 * i ) + 2] - dphi[i][qp]( 2 ) );
            Bl_np1[2][( 3 * i ) + 2] =
                dphi[i][qp]( 2 ) + ( one3 * ( avgTraceTerm[( 3 * i ) + 2] - dphi[i][qp]( 2 ) ) );
            Bl_np1[3][( 3 * i ) + 1] = dphi[i][qp]( 2 );
            Bl_np1[3][( 3 * i ) + 2] = dphi[i][qp]( 1 );
            Bl_np1[4][( 3 * i ) + 0] = dphi[i][qp]( 2 );
            Bl_np1[4][( 3 * i ) + 2] = dphi[i][qp]( 0 );
            Bl_np1[5][( 3 * i ) + 0] = dphi[i][qp]( 1 );
            Bl_np1[5][( 3 * i ) + 1] = dphi[i][qp]( 0 );
        }

        for ( int i = 0; i < 6; i++ )
            for ( unsigned int j            = 0; j < ( 3 * num_nodes ); j++ )
                materialStiffnessTemp[i][j] = 0.0;

        for ( unsigned int i = 0; i < ( 3 * num_nodes ); i++ )
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) materialStiffness[i][j] = 0.0;

        for ( int i = 0; i < 6; i++ )
            for ( int j = 0; j < 6; j++ ) materialMatrix[i][j] = 0.0;

        d_materialModel->getConstitutiveMatrix( constitutiveMatrix );

        for ( int i = 0; i < 6; i++ ) {
            for ( int j = 0; j < 6; j++ ) {
                materialMatrix[i][j] = constitutiveMatrix[( 6 * i ) + j];
            }
        }

        for ( int i = 0; i < 6; i++ ) {
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                for ( int k = 0; k < 6; k++ ) {
                    materialStiffnessTemp[i][j] += ( materialMatrix[i][k] * Bl_np1[k][j] );
                }
            }
        }

        for ( unsigned int i = 0; i < ( 3 * num_nodes ); i++ ) {
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                for ( int k = 0; k < 6; k++ ) {
                    materialStiffness[i][j] += ( Bl_np1[k][i] * materialStiffnessTemp[k][j] );
                }
            }
        }

        for ( unsigned int j = 0; j < num_nodes; j++ ) {
            for ( int d1 = 0; d1 < 3; d1++ ) {

                /*          // Strain Displacement Relation For Columns (Voigt Notation)
                          double uStrain[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

                          if(d1 == 0) {
                            uStrain[0] = dphi[j][qp](0); //xx
                            uStrain[5] = (0.5*dphi[j][qp](1)); //xy
                            uStrain[4] = (0.5*dphi[j][qp](2)); //xz
                          }

                          if(d1 == 1) {
                            uStrain[5] = (0.5*dphi[j][qp](0)); //yx
                            uStrain[1] = dphi[j][qp](1); //yy
                            uStrain[3] = (0.5*dphi[j][qp](2)); //yz
                          }

                          if(d1 == 2) {
                            uStrain[4] = (0.5*dphi[j][qp](0)); //zx
                            uStrain[3] = (0.5*dphi[j][qp](1)); //zy
                            uStrain[2] = dphi[j][qp](2); //zz
                          }

                          for(int c = 0; c < 3; c++) {
                            uStrain[c] = uStrain[c] + ((avgTraceTerm[(3*j) + d1] -
                   dphi[j][qp](d1))/3.0);
                          }//end for c

                          // Stress-Strain Relation in Voigt Notation
                          double uStress[6];

                          for(int r = 0; r < 6; r++)  {
                            uStress[r] = 0;
                            for(int c = 0; c < 6; c++) {
                              uStress[r] += (constitutiveMatrix[(6*r) + c]*uStrain[c]);
                            }//end for c
                          }//end for r

                          // Convert Voigt Notation to Matrix Notation
                          double T[3][3];

                          T[0][0] = uStress[0]; //xx
                          T[0][1] = uStress[5]; //xy
                          T[0][2] = uStress[4]; //xz

                          T[1][0] = uStress[5]; //yx
                          T[1][1] = uStress[1]; //yy
                          T[1][2] = uStress[3]; //yz

                          T[2][0] = uStress[4]; //zx
                          T[2][1] = uStress[3]; //zy
                          T[2][2] = uStress[2]; //zz
                */
                for ( unsigned int k = 0; k < num_nodes; k++ ) {
                    for ( int d2 = 0; d2 < 3; d2++ ) {

                        /*              // Strain Displacement Relation For Rows (Voigt Notation)
                                      double wStrain[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

                                      if(d2 == 0) {
                                        wStrain[0] = dphi[k][qp](0); //xx
                                        wStrain[5] = (0.5*dphi[k][qp](1)); //xy
                                        wStrain[4] = (0.5*dphi[k][qp](2)); //xz
                                      }

                                      if(d2 == 1) {
                                        wStrain[5] = (0.5*dphi[k][qp](0)); //yx
                                        wStrain[1] = dphi[k][qp](1); //yy
                                        wStrain[3] = (0.5*dphi[k][qp](2)); //yz
                                      }

                                      if(d2 == 2) {
                                        wStrain[4] = (0.5*dphi[k][qp](0)); //zx
                                        wStrain[3] = (0.5*dphi[k][qp](1)); //zy
                                        wStrain[2] = dphi[k][qp](2); //zz
                                      }

                                      for(int c = 0; c < 3; c++) {
                                        wStrain[c] = wStrain[c] + ((avgTraceTerm[(3*k) + d2] -
                           dphi[k][qp](d2))/3.0);
                                      }//end for c

                                      // Convert Voigt Notation to Matrix Notation
                                      double Sw[3][3];

                                      Sw[0][0] = wStrain[0]; //xx
                                      Sw[0][1] = wStrain[5]; //xy
                                      Sw[0][2] = wStrain[4]; //xz

                                      Sw[1][0] = wStrain[5]; //yx
                                      Sw[1][1] = wStrain[1]; //yy
                                      Sw[1][2] = wStrain[3]; //yz

                                      Sw[2][0] = wStrain[4]; //zx
                                      Sw[2][1] = wStrain[3]; //zy
                                      Sw[2][2] = wStrain[2]; //zz

                                      // trace{TS}
                                      double tmp = 0;
                                      for(int r = 0; r < 3; r++) {
                                        for(int c = 0; c < 3; c++) {
                                          tmp += (T[r][c]*Sw[c][r]);
                                        }
                                      }
                        */
                        // elementStiffnessMatrix[(3*k) + d2][(3*j) + d1] += (JxW[qp]*tmp);
                        elementStiffnessMatrix[( 3 * k ) + d2][( 3 * j ) + d1] +=
                            ( JxW[qp] * materialStiffness[( 3 * k ) + d2][( 3 * j ) + d1] );

                    } // end for d2
                }     // end for k
            }         // end for d1
        }             // end for j

        d_materialModel->postLinearGaussPointOperation();
    } // end for qp

    d_materialModel->postLinearElementOperation();
}

void MechanicsLinearElement::apply_Normal()
{
    const std::vector<Real> &JxW = ( *d_JxW );

    const std::vector<std::vector<RealGradient>> &dphi = ( *d_dphi );

    std::vector<std::vector<double>> &elementStiffnessMatrix = ( *d_elementStiffnessMatrix );

    d_fe->reinit( d_elem );

    d_materialModel->preLinearElementOperation();

    const unsigned int num_nodes = d_elem->n_nodes();

    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
        d_materialModel->preLinearGaussPointOperation();

        double *constitutiveMatrix;
        double Bl_np1[6][24], materialMatrix[6][6];
        double materialStiffness[24][24], materialStiffnessTemp[6][24];

        for ( int i = 0; i < 6; i++ ) {
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                Bl_np1[i][j] = 0.0;
            }
        }

        for ( unsigned int i = 0; i < num_nodes; i++ ) {
            Bl_np1[0][( 3 * i ) + 0] = dphi[i][qp]( 0 );
            Bl_np1[1][( 3 * i ) + 1] = dphi[i][qp]( 1 );
            Bl_np1[2][( 3 * i ) + 2] = dphi[i][qp]( 2 );
            Bl_np1[3][( 3 * i ) + 1] = dphi[i][qp]( 2 );
            Bl_np1[3][( 3 * i ) + 2] = dphi[i][qp]( 1 );
            Bl_np1[4][( 3 * i ) + 0] = dphi[i][qp]( 2 );
            Bl_np1[4][( 3 * i ) + 2] = dphi[i][qp]( 0 );
            Bl_np1[5][( 3 * i ) + 0] = dphi[i][qp]( 1 );
            Bl_np1[5][( 3 * i ) + 1] = dphi[i][qp]( 0 );
        }

        for ( int i = 0; i < 6; i++ )
            for ( unsigned int j            = 0; j < ( 3 * num_nodes ); j++ )
                materialStiffnessTemp[i][j] = 0.0;

        for ( unsigned int i = 0; i < ( 3 * num_nodes ); i++ )
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) materialStiffness[i][j] = 0.0;

        for ( int i = 0; i < 6; i++ )
            for ( int j = 0; j < 6; j++ ) materialMatrix[i][j] = 0.0;

        d_materialModel->getConstitutiveMatrix( constitutiveMatrix );

        for ( int i = 0; i < 6; i++ ) {
            for ( int j = 0; j < 6; j++ ) {
                materialMatrix[i][j] = constitutiveMatrix[( 6 * i ) + j];
            }
        }

        for ( int i = 0; i < 6; i++ ) {
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                for ( int k = 0; k < 6; k++ ) {
                    materialStiffnessTemp[i][j] += ( materialMatrix[i][k] * Bl_np1[k][j] );
                }
            }
        }

        for ( unsigned int i = 0; i < ( 3 * num_nodes ); i++ ) {
            for ( unsigned int j = 0; j < ( 3 * num_nodes ); j++ ) {
                for ( int k = 0; k < 6; k++ ) {
                    materialStiffness[i][j] += ( Bl_np1[k][i] * materialStiffnessTemp[k][j] );
                }
            }
        }

        for ( unsigned int j = 0; j < num_nodes; j++ ) {
            for ( int d1 = 0; d1 < 3; d1++ ) {

                /*          // Strain Displacement Relation For Columns (Voigt Notation)
                          double uStrain[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

                          if(d1 == 0) {
                            uStrain[0] = dphi[j][qp](0); //xx
                            uStrain[5] = (0.5*dphi[j][qp](1)); //xy
                            uStrain[4] = (0.5*dphi[j][qp](2)); //xz
                          }

                          if(d1 == 1) {
                            uStrain[5] = (0.5*dphi[j][qp](0)); //yx
                            uStrain[1] = dphi[j][qp](1); //yy
                            uStrain[3] = (0.5*dphi[j][qp](2)); //yz
                          }

                          if(d1 == 2) {
                            uStrain[4] = (0.5*dphi[j][qp](0)); //zx
                            uStrain[3] = (0.5*dphi[j][qp](1)); //zy
                            uStrain[2] = dphi[j][qp](2); //zz
                          }

                          // Stress-Strain Relation in Voigt Notation
                          double uStress[6];

                          for(int r = 0; r < 6; r++)  {
                            uStress[r] = 0;
                            for(int c = 0; c < 6; c++) {
                              uStress[r] += (constitutiveMatrix[(6*r) + c]*uStrain[c]);
                            }//end for c
                          }//end for r

                          // Convert Voigt Notation to Matrix Notation
                          double T[3][3];

                          T[0][0] = uStress[0]; //xx
                          T[0][1] = uStress[5]; //xy
                          T[0][2] = uStress[4]; //xz

                          T[1][0] = uStress[5]; //yx
                          T[1][1] = uStress[1]; //yy
                          T[1][2] = uStress[3]; //yz

                          T[2][0] = uStress[4]; //zx
                          T[2][1] = uStress[3]; //zy
                          T[2][2] = uStress[2]; //zz
                */
                for ( unsigned int k = 0; k < num_nodes; k++ ) {
                    for ( int d2 = 0; d2 < 3; d2++ ) {

                        /*              // Strain Displacement Relation For Rows (Voigt Notation)
                                      double wStrain[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

                                      if(d2 == 0) {
                                        wStrain[0] = dphi[k][qp](0); //xx
                                        wStrain[5] = (0.5*dphi[k][qp](1)); //xy
                                        wStrain[4] = (0.5*dphi[k][qp](2)); //xz
                                      }

                                      if(d2 == 1) {
                                        wStrain[5] = (0.5*dphi[k][qp](0)); //yx
                                        wStrain[1] = dphi[k][qp](1); //yy
                                        wStrain[3] = (0.5*dphi[k][qp](2)); //yz
                                      }

                                      if(d2 == 2) {
                                        wStrain[4] = (0.5*dphi[k][qp](0)); //zx
                                        wStrain[3] = (0.5*dphi[k][qp](1)); //zy
                                        wStrain[2] = dphi[k][qp](2); //zz
                                      }

                                      // Convert Voigt Notation to Matrix Notation
                                      double Sw[3][3];

                                      Sw[0][0] = wStrain[0]; //xx
                                      Sw[0][1] = wStrain[5]; //xy
                                      Sw[0][2] = wStrain[4]; //xz

                                      Sw[1][0] = wStrain[5]; //yx
                                      Sw[1][1] = wStrain[1]; //yy
                                      Sw[1][2] = wStrain[3]; //yz

                                      Sw[2][0] = wStrain[4]; //zx
                                      Sw[2][1] = wStrain[3]; //zy
                                      Sw[2][2] = wStrain[2]; //zz

                                      // trace{TS}
                                      double tmp = 0;
                                      for(int r = 0; r < 3; r++) {
                                        for(int c = 0; c < 3; c++) {
                                          tmp += (T[r][c]*Sw[c][r]);
                                        }
                                      }
                        */
                        // elementStiffnessMatrix[(3*k) + d2][(3*j) + d1] += (JxW[qp]*tmp);
                        elementStiffnessMatrix[( 3 * k ) + d2][( 3 * j ) + d1] +=
                            ( JxW[qp] * materialStiffness[( 3 * k ) + d2][( 3 * j ) + d1] );

                    } // end for d2
                }     // end for k
            }         // end for d1
        }             // end for j

        d_materialModel->postLinearGaussPointOperation();
    } // end for qp

    d_materialModel->postLinearElementOperation();
}
}
}
