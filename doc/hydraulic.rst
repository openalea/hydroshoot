===================
Hydraulic structure
===================

The *hydraulic* module computes the distribution of xylem potential across plant segments by analogy to Ohm’s law
(:numref:`_fig_1_hydraulic`)

.. _fig_1_hydraulic:

.. figure:: fig/fig_1_hydraulic.png

    Illustration of the variables required to calculate the hydraulic structure: xylem flux, (:math:`F`),
    conductivity to xylem flux (:math:`K`), xylem pressure at upper (downstream) and lower (upstream) extremities
    of the conducting element (respectively :math:`H_u` and :math:`H_l`), and the length of the segment (:math:`L`).

Thereby, xylem flux (:math:`F \ [kg \ s^{-1}]`) across the hydraulic segment of length (:math:`L \ [m]`)
is driven by the difference of xylem pressures across this segment (:math:`H_u-H_l \ [MPa]`)
and regulated by segment's conductivity to xylem flow (:math:`K \ [kgs^{-1} \ m \ MPa^{-1}]`):

.. math::
    F = - K \dot \frac{H_u - H_l}{L}

Xylem conductivity varies with water potential as a result of xylem cavitation under water deficit
**(Tyree and Sperry, 1989)**. This relationship is described in HydroShoot as:

.. math::
    K = K_{max} \ \frac{1} {\left( 1 + \left( \frac{\Psi}{\Psi_{crit, \ stem}} \right) ^{c_{x1}} \right)}

where
:math:`K_{max} \ [kg \ s^{-1} \ m \ MPa^{-1}]` is the maximum conductivity of the segment,
:math:`\Psi \ [MPa]` is the arithmetic mean of xylem potential of the segment,
:math:`\Psi_{crit, \ stem} \ [MPa]` is a shape parameter defining **XXXXXX**, and
:math:`c_{x1} \ [-]` is shape parameter defining the steepness of the response of :math:`K` to :math:`\Psi`.

Finally, :math:`K_{max}` is estimated empirically as proposed by **Tyree and Zimmermann (2002)** as:

..math::
    K_{max} = c_{x2} \dot D^{c_{x3}}

where
:math:`D \ [m]` is the average diameter of the segment, and
:math:`c_{x2}` and
:math:`c_{x3}` are shape parameters, mostly given within the ranges of [2.5, 2.8] and [2.0, 5.0], respectively.

The last two equations apply to all conducting segments (not leaves blades). Xylem potential of the upper extremity
of the petiole is assumed equal to that of the lumped leaf water potential :math:`\Psi_{leaf} \ [MPa]`.

References
----------
Tyree M, Sperry J. 1989.
    Vulnerability of xylem to cavitation and embolism.
    Annual review of plant physiology and plant molecular biology 40: 19–38.
Tyree M, Zimmermann M. 2002.
    Xylem structure and the ascent of sap, Springer Series in Wood Science.


