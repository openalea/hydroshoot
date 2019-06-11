============
Gas exchange
============

The exchange module computes the rates of net carbon assimilation and transpiration per unit leaf surface area
(respectively :math:`A_n` and :math:`E`) as a function of micrometeorological conditions and leaf water status.
The calculations use the analytical solution proposed by **Yin et al. (2009)** for coupling the C_3 photosynthesis
model of **Farquhar et al. (1980)** to the stomatal conductance model of **Ball et al. (1987)**. This coupling
allows stomatal conductance (to both CO2 and water vapor) to respond to environmental stimuli
(temperature and irradiance) via photosynthesis. It is based on Fick’s first law of diffusion, whereby the net
assimilation :math:`A_n`, the stomatal conductance to CO2 :math:`g_{s, \ CO_2}`, and the mesophyll conductance
:math:`g_m` are implied. The solution is based on the three following equations (following **Evers et al. 2010**
supporting information):

.. math::
    A_n = \frac{(C_c - \Gamma) \dot x_1}{C_c + x_2} - R_d
    C_c = C_i - \frac{A_n}{g_m}
    g_{s, \ CO_2} = g_{s0, \ CO_2} + m_0 \dot \frac{A_n + R_d}{C_i - \Gamma} \dot f_w
    g_{s, \ CO_2} = \frac{A_n}{C_a - C_i - A_n \dot r_{tb}}

where
:math:`A_n \ [\mu mol \ m^{-2} \ s^{-1}]` is net carbon assimilation rate,
:math:`R_d \ [\mu mol \ m^{-2} \ s^{-1}]` is mitochondrial respiration in the light,
:math:`\Gamma \ [\mu bar]` is CO_2 compensation point in the absence of mitochondrial respiration,
:math:`x_1 \ [\mu mol \ m^{-2} \ s^{-1}]` and :math:`x_2 \ [\mu bar]` are intermediate parameters,
:math:`g_m \ [\mu mol \ m^{-2} \ s^{-1} \ {\mu bar}^{-1}]` is mesophyll conductance for CO_2 diffusion,
:math:`g_{s, \ CO_2} \ [mol \ m^{-2} \ s^{-1} \ {\mu bar}^{-1}]` is stomatal conductance to CO_2 diffusion,
:math:`g_{s0, \ CO_2} \ [mol \ m^{-2} \ s^{-1} \ {\mu bar}^{-1}]` is the residual stomatal conductance to CO_2 diffusion,
:math:`m_0 \ [-]` is a shape parameter regulating the slope between :math:`A_n` and :math:`g_{s, \ CO_2}`,
:math:`f_w \ [-]` is a dimensionless function representing the response of stomatal conductance to soil or plant water status,
:math:`r_{tb} \ [m^2 \ s \ \mu bar \ {\mu mol}^{-1})]` is the combined turbulence and boundary layer resistance for CO2,
:math:`C_a \ [\mu bar]` is air CO_2 partial pressure,
:math:`C_i \ [\mu bar]` is intercellular CO_2 partial pressure, and
:math:`C_c \ [\mu bar]` is chloroplast CO_2 partial pressure.


Net carbon assimilation rate per unit leaf area
===============================================

Variable intra-canopy photosynthetic capacities
-----------------------------------------------

Transpiration rate per unit leaf area
=====================================

The transpiration rate :math:`E \ [mol \ m^{-2} s^(-1)]` is calculated as:

.. math::
    E = \frac{1}{\frac{1}{g_{b, \ H_2O} + \frac{1}{1.6 \dot g_{s, \ CO_2}}} left( \frac{VPD}{P_a} right)


where
:math:`P_a \ [kPa]` is the atmospheric pressure and
:math:`g_{b, \ H_2O} \ [mol \ m^{-2} s^(-1)]` is the boundary layer conductance to water vapor derived from
**Nobel (2005)** as:

.. math::
    g_{b, \ H_2O} = \frac{D_{H_2O} \dot P_v}{R \dot T_{leaf} \dot \Delta x}

with

.. math::
    D_{H_2O} = D_(H_2O, 0) \frac{P_a}{P_v} left( \frac{T_{leaf}}{273} right)^1.8

where
:math:`D_{H_2O}` is the diffusion coefficient of H2O in the air at 0 :math:`^\circ C` (:math:`2.13e-5 \ m^2 s^{-1}`),
:math:`P_a \ [MPa]` is the ambient air pressure at 0 :math:`^\circ C`,
:math:`P_v \ [MPa]` is water vapor partial pressure, and
:math:`\Delta x \ [m]` is the thickness of the boundary layer defined as **(Nobel 2005)**:

.. math::
    \Delta x = 0.004 \sqrt{\frac{l}{v}}

where
:math:`l \ [m]` is the mean length of the leaf in the downwind direction (set to 70% of blade length), and
:math:`v \ [m \ s^{-1}]` is wind speed in the vicinity of the leaf.


The impact of water stress on stomatal conductance (i.e. via the :math:`f_w` function) is calculated using one of the
following options:

.. math::
    f_w =   \left \{
                \begin{array}{11}
                    \frac{1}{1+(VPD/D_0)}                                    &   (a) \\
                    \frac{1}{1+(\frac{\Psi_{leaf}}{\Psi_{crit, \ leaf}})^n}   &   (b) \\
                    \frac{1}{1+(\frac{\Psi_{soil}}{\Psi_{crit, \ leaf}})^n}   &   (c) \\
                \end{array}
            \right.

where
:math:`VPD \ [kPa]` is vapor pressure deficit (between the leaf and the air),
:math:`D_0 \ [kPa]` shape parameter,
:math:`Psi_{leaf} \ [MPa]` leaf bulk xylem potential,
:math:`Psi_{soil} \ [MPa]` soil bulk water potential (assumed equal to xylem potential at the base of the shoot), and
:math:`Psi_{crit, leaf} \ [MPa]` leaf water potential at which stomatal conductance reduces to half its maximum value.


In case the option :math:`a` is used, stomatal conductance reduction is considered independent from the soil water
status (i.e. following **Leuning, 1995**). In contrast, Both options :math:`b` and `c` allows simulating stomatal
conductance as a function of the leaf xylem potential (i.e. regarding shoot hydraulic structure) or the soil
water potential (i.e. disregarding the hydraulic structure of the shoot), respectively.


As Farquhar’s model has been thoroughly detailed in literature, its description is given in Appendix I.





g_(s,CO_2 ) is calculated according to **Yin et al. (2009)** as:

g_(s,CO_2 )=g_(s0,CO_2 )+

where g_(s0,CO_2 ) is the residual stomatal conductance to 〖CO〗_2 [〖mol〗_(〖CO〗_2 )  m^(-2)  s^(-1) ],

  mo〖l_(CO_2 )〗^(-1) ], m_0 is a dimensionless shape parameter, and f_w is a dimensionless function representing the response of g_(s,CO_2 ) to air water vapor deficit (VPD,kPa). f_w is deduced from the stomatal conductance model of Leuning (1995) as:

f_w=1/((1+VPD/D_0 ) )	(Eq. 5a)

where D_0 is a scaling parameter [kPa].

Eq. 5a does not account for stomatal sensitivity to soil water deficit (“remote” approach) or local leaf water potential (“local” approach). Tuzet et al. (2003) and Leuning et al. (2004) suggested to express f_w as a function of the local Ψ_leaf. This function is implemented in HydroShoot following Nikolov (1995):

f_w=1/((1+(Ψ_leaf/Ψ_(crit,leaf) )^n ) )	(Eq. 5b)

where Ψ_(crit,leaf) is a critical leaf water potential threshold [MPa] at which stomatal conductance is reduced by 50%, and n is a shape parameter [-]. The same last equation is used to express the dependency of g_(s0,CO_2 ) on the remote soil water potential (Ψ_soil):

f_w=1/((1+(Ψ_soil/Ψ_(crit,leaf) )^n ) )	(Eq. 5c)


Intra-canopy variability in photosynthetic capacity

Leaf photosynthetic traits (maximum carboxylation rate V_cmax , maximum electron transport rate J_max, triose-phosphate transport rate TPU and R_d; cf. Appendix I) have been shown to strongly vary within the plant canopy so that to increase light-saturated net assimilation rate with increasing solar irradiance availability throughout the canopy (Niinemets et al., 2014). HydroShoot accounts for this variability by considering leaf nitrogen content per unit leaf surface area (N_a,g_N  m^(-2)) as the pivotal trait to determine the photosynthetic capacity of leaves (Prieto et al., 2012) as follows:

P^25=S_(N_a ) N_a-b_(N_a )	(Eq. 10)

where P^25 is the value at 25 °C for any of the rates V_cmax, J_max, TPU and R_d (given as inputs), and S_(N_a ) [μmol_(CO_2 )  〖g_N〗^(-1)  s^(-1) ] and b_(N_a ) [μmol_(CO_2 )  m^(-2)  s^(-1) ] are the slope and the intercept of the linear relationship with N_a specific to each rate. N_a is calculated as the product of nitrogen content per unit leaf dry mass N_m [g_N  〖g_drymatter〗^(-1) ] and leaf dry mass per area LMA [g_drymatter  m^(-2) ]. N_m linearly varies with plant age, expressed as the thermal time cumulated from budburst (input of the model), and LMA is determined by leaf exposure to light during the last past days (Prieto et al., 2012). This is expressed respectively in the two following equations:

N_m=a_N ∑_(i=budburst)^d▒(max(0,T_(air,i)-T_b )) +b_N	(Eq. 11)
LMA=a_M ln(PPFD_10 )+b_M	(Eq. 12)

where T_(air,i) is the mean temperature of the day i [°C] and T_b is the base temperature (minimum required for growth) [°C], set to 10°C for grapevine and used for the calculation of thermal time since budburst, a_N [g_N  〖g_drymatter〗^(-1)°Cd^(-1) ] and b_N [g_N  〖g_drymatter〗^(-1) ] are the slope and intercept of the linear relationship between N_m and accumulated thermal time since budburst, PPFD_10 [mol_photon  m^(-2)  d^(-1) ] is the cumulative photosynthetic photon flux density irradiance intercepted by the leaf (output of the energy module) averaged over the past 10 days, a_M [g_drymatter  m^(-2) ] and b_M [g_drymatter  m^(-2) ] are the slope and intercept of the linear relationship between LMA and the logarithm of PPFD_10.

Finally, this module was provided with a photoinhibition model as this phenomenon is frequently reported to affect grapevines under combined heat and water stresses (Correia et al., 1990; Flexas and Medrano, 2002; Lovisolo et al., 2010). The simple photoinhibition model implemented in HydroShoot is detailed in Appendix II and assumes that combined heat and water stresses inhibit photosynthesis by reducing the electron transport rate (cf. J in Eq. A6) as the result of an increase of deactivation energy ΔH_d (cf. equations A9 and A10).


yin et al
farquhar
Evers et al. (2010)
leuning 1995
nobel 2005