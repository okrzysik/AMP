#ifndef included_OxideModel
#define included_OxideModel


#include <string>


namespace AMP{
namespace TimeIntegrator{


/*!
  @brief This class provides routines for oxide calculations at a point
 */
class OxideModel
{
public:

    /*!
     * @brief Function to return the equiblibrium concentration at the phase 
     *  boundaries in zirconium.
     *
     * @param T     Temperature to use for the equiblibrium (K)
     * @param C     Returned equilibrium concentrations at the phase boundaries.
     *              The output will be a 1x6 array:
     *              C[0] - The concentration in the oxide at the oxide-gas interface (g/cm^3)
     *              C[1] - The concentration in the oxide at the oxide-alpha interface (g/cm^3)
     *              C[2] - The concentration in the alpha at the alpha-oxide interface (g/cm^3)
     *              C[3] - The concentration in the alpha at the alpha-beta interface (g/cm^3)
     *              C[4] - The concentration in the beta at the beta-alpha interface (g/cm^3)
     *              C[5] - The concentration in the beta at the beta-bulk interface (g/cm^3)
     */
    static void get_equilibrium_concetration( double T, double *C );


    /*!
     * @brief Function to return the diffusion coefficients in the layers
     *
     * @param T     Temperature to use for the diffusion coefficients (K)
     * @param D     Returned diffusion coefficients in the layers (cm^2/s).
     *              The output will be a 1x3 array:
     *              D[0] - The diffusion coefficient in the oxide layer
     *              D[1] - The diffusion coefficient in the alpha layer
     *              D[2] - The diffusion coefficient in the beta layer
     */
    static void get_diffusion_coefficients( double T, double *D );


    /*!
     * @brief Function to integrate the oxide over a fixed timestep.
     *  This routine will take the current solution and use an internal
     *  time stepping to advance the solution to the next point in time.
     * @return value is the number of internal timesteps taken.
     * @param[in] dT    Timestep (s)
     * @param[in] N     Number of layers
     * @param[in] N2    Number of points in each layer
     * @param[in] x0    Initial boundary position (cm) (N+1)
     * @param[in] Cb    Concentration at boundaries (g/cm^3) (2xN)
     * @param[in] C0    Initial concentration of each layer (g/cm^3) {N}(1xN2)
     * @param[in] D     Diffusion coefficient of each layer (cm^2/s) (N)
     * @param[out] C1   Final concentration of each layer (g/cm^3) {N}(1xN2)
     * @param[out] x1   Final boundary position (cm) (N+1)
     * @param[out] v1   Final boundary velocity (cm/s) (N+1)
     */
    static int integrateOxide( double dT, int N, const int *N2, const double *x0, const double *Cb, 
        double **C0, const double *D, double **C1, double *x1, double *v1 );


private:

    /*!
     * @brief This function computes the velocity of the bounary layers.
     * @param[in] N     Number of layers 
     * @param[in] x     Boundary position (N+1)
     * @param[in] Cb    Concentration at boundaries (2*N)
     * @param[in] Nl    Number of zones in each layer (N)
     * @param[in] C     Concentration of each layer ( N x (Nl(i)) )
     * @param[in] D     Diffusion coefficient of each layer (N)
     * @param[out] v    Velocity of each boundary interface (N+1)
     */
    static void computeVelocity( const int N, const double *x, const double *Cb, 
        const int *Nl, double const* const* C, const double *D, double *v );


    /* This function computes the maximum timestep based on diffusion
     * @param[in] N     Number of zones
     * @param[in] x     Boundary position
     * @param[in] Cb    Concentration at boundaries
     * @param[in] C     Concentration of each layer (N)
     * @param[in] D     Diffusion coefficient
     */
    static double computeDiffustionTimestep( const int N, const double x[2], 
        const double Cb[2], const double *C, const double D );


    /* This function performs a linear solve to get the future oxygen concentration
     * given the current oxygen concentration, the current and future boundary 
     * velocities, and the current boundary positions.
     * @param[in] N     Number of zones
     * @param[in] dt    Timestep
     * @param[in] x0    Initial boundary position
     * @param[in] v0    Velocity of boundaries
     * @param[in] Cb    Concentration at boundaries
     * @param[in] C0    Initial concentration (N)
     * @param[in] D     Diffusion coefficient
     * @param[out] C1   Final concentration (N)
     * @param[out] x1   Final boundary position
     */
    static void solveLinearDiffusionLayer( const int N, const double dt, const double x0[2], 
        const double v[2], const double Cb[2], const double *C0,  const double D,
        double *C1, double *x1 );


};


}
}

#endif
