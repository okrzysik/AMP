#ifndef included_AMP_Constants
#define included_AMP_Constants

#include "AMP/utils/Units.h"


namespace AMP::Constants {


//! pi
constexpr double pi = 3.1415926535897932;

//! Newtonian constant of gravitation
constexpr Units G( "m^3 kg^-1 s^-2", 6.6743015e-11 );

//! speed of light in vacuum
constexpr Units c( "m*s^-1", 299792458 );

//! elementary charge
constexpr Units e( "C", 1.602176634e-19 );

//! Planck constant
constexpr Units h( "J*Hz^-1", 6.62607015e-34 );

//! reduced Planck constant: h/(2*pi)
constexpr Units h_bar( "J*s", 1.054571817e-34 );

//! Coulomb constant
constexpr Units ke( "N*m^2*C^-2", 8.987551792314e9 );

//! Boltzmann constant
constexpr Units kB( "J*K^-1", 1.380649e-23 );

//! electron mass
constexpr Units me( "kg", 9.109383701528e-31 );

//! proton mass
constexpr Units mp( "kg", 1.6726219236951e-27 );

//! neutron mass
constexpr Units mn( "kg", 1.6749274980495e-27 );

//! muon mass
constexpr Units m_mu( "kg", 1.88353162742e-28 );

//! tau mass
constexpr Units m_tau( "kg", 3.1675421e-27 );

//! top quark mass
constexpr Units m_t( "kg", 3.078453e-25 );

//! vacuum electric permittivity
constexpr Units epsilon0( "F*m^-1", 8.854187812813e-12 );

//! vacuum magnetic permeability
constexpr Units mu0( "N*A^-2", 1.2566370621219e-6 );

//! characteristic impedance of vacuum: mu0*c
constexpr Units Z0( "ohm", 376.73031366857 );

//! Stefan–Boltzmann constant: pi^2*kB^4/(60*k_bar^3*c^2)
constexpr Units sigma( "W*m^-2*K^-4", 5.670374419e-8 );

//! first radiation constant: 2*pi*h*c^2
constexpr Units ca( "W*m^2", 3.741771852e-16 );

//! first radiation constant for spectral radiance: 2*h*c/sr
constexpr Units c1L( "W*m^2*sr^-1", 1.1910429723971884140794892e-16 );

//! second radiation constant: h*c/kB
constexpr Units c2( "m*K", 1.438776877e-2 );

//! Wien wavelength displacement law constant
constexpr Units b( "m*K", 2.897771955e-3 );

//! Wien frequency displacement law constant
constexpr Units bf( "Hz*K^-1", 5.878925757e10 );

//! Wien entropy displacement law constant
constexpr Units b_entropy( "m*K", 3.002916077e-3 );

//! conductance quantum: 2*e^2/h
constexpr Units G0( "S", 7.748091729e-5 );

//! inverse conductance quantum: h/(2*e^2)
constexpr Units G0_1( "ohm", 12906.40372 );

//! von Klitzing constant: h/e^2
constexpr Units RK( "ohm", 25812.80745 );

//! Josephson constant: 2*e/h
constexpr Units HJ( "Hz*V^-1", 483597.8484e9 );

//! magnetic flux quantum: h/(2*e)
constexpr Units I0( "Wb", 2.067833848e-15 );

//! fine-structure constant: e^2/(4*pi*eplison0*h_bar*c)
constexpr Units alpha( "", 7.297352569311e-3 );

//! inverse fine-structure constant: (4*pi*eplison0*h_bar*c)/e^2
constexpr Units alpha_1( "", 137.03599908421 );

//! proton-to-electron mass ratio: mp/me
constexpr Units mp_me( "", 1836.1526734311 );

//! W-to-Z mass ratio: m_W/m_Z
constexpr Units mW_mZ( "", 0.8815317 );

//! weak mixing angle: sin^2(theta)=1-(m_W/m_Z)^2
constexpr Units theta_W( "", 0.2229030 );

//! electron g-factor
constexpr Units g_e( "", -2.0023193043625635 );

//! muon g-factor
constexpr Units g_mu( "", -2.002331841813 );

//! proton g-factor
constexpr Units g_p( "", 5.585694689316 );

//! quantum of circulation: h_bar/(2*me)
constexpr Units quantum_circulation( "m^2*s^-1", 3.636947551611e-4 );

//! Bohr magneton: e*h_bar/(2*me)
constexpr Units mu_B( "J*T^-1", 9.274010078328e-24 );

//! nuclear magneton: e*h_bar/(2*mp)
constexpr Units mu_N( "J*T^-1", 5.050783746115e-27 );

//! classical electron radius: e^2*k_e/(me*c^2)
constexpr Units r_e( "m", 2.817940326213e-15 );

//! Thomson cross section: (8*pi/3)*r_e^2
constexpr Units sigma_e( "m^2", 6.652458732160e-29 );

//! Bohr radius: h_bar^2/(k_e*me*e^2) or r_e/alpha^2
constexpr Units a0( "m", 5.2917721090380e-11 );

//! Hartree energy: alpha^2*c^2*me
constexpr Units Eh( "J", 4.359744722207185e-18 );

//! Rydberg unit of energy: Eh/2
constexpr Units Ry( "J", 2.179872361103542e-18 );

//! Rydberg constant: alpha^2*me*c/(2*h)
constexpr Units R_inf( "m^-1", 10973731.56816021 );

//! Fermi coupling constant: Gf/(h_bar*c^3)
constexpr Units Fermi_coupling( "GeV^-2", 1.16637876e-5 );

//! Avogadro constant
constexpr Units NA( "mol^-1", 6.02214076e23 );

//! molar gas constant: NA*kB
constexpr Units R( "J*mol^-1*K^-1", 8.31446261815324 );

//! Faraday constant: NA*e
constexpr Units F( "C*mol^-1", 96485.3321233100184 );

//! molar Planck constant: NA*h
constexpr Units molar_Plank( "J*s*mol^-1", 3.9903127128934314e-10 );

//! atomic mass of carbon-12
constexpr Units m_12C( "kg", 1.9926468799260e-26 );

//! molar mass of carbon-12
constexpr Units M_12C( "kg*mol^-1", 11.999999995836e-3 );

//! atomic mass constant: m_12C/12
constexpr Units m_u( "kg", 1.6605390666050e-27 );

//! molar mass constant: M_12C/12
constexpr Units M_u( "kg*mol^-1", 0.9999999996530e-3 );

//! molar volume of silicon
constexpr Units Vm_Si( "m^3*mol^-1", 1.20588319960e-5 );

//! hyperfine transition frequency of 133Cs
constexpr Units dnu_Cs( "Hz", 9192631770 );


} // namespace AMP::Constants

#endif
