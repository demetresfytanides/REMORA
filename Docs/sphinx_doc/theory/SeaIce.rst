.. role:: cpp(code)
   :language: c++

.. _sec:SeaIce-equations:

Sea-Ice Equations
=================

What the model solves for
-------------------------

The unknowns are two-dimensional fields, one value per ocean column:

+-----------------------------------------------+---------------------------------------------------+-----------------------------+
| Symbol                                        | Unknown                                           | Solved by                   |
+===============================================+===================================================+=============================+
| :math:`A`                                     | Ice concentration, ``Aice``                       | Thermodynamics, advection   |
+-----------------------------------------------+---------------------------------------------------+-----------------------------+
| :math:`h`                                     | Ice thickness, cell average [m]                   | Thermodynamics, advection   |
+-----------------------------------------------+---------------------------------------------------+-----------------------------+
| :math:`h_s`                                   | Snow thickness, cell average [m]                  | Thermodynamics, advection   |
+-----------------------------------------------+---------------------------------------------------+-----------------------------+
| :math:`h_m`                                   | Melt-pond thickness, cell average [m]             | Thermodynamics, advection   |
+-----------------------------------------------+---------------------------------------------------+-----------------------------+
| :math:`T_i`                                   | Interior ice temperature, ``Tice`` [C]            | Thermodynamics, advection   |
+-----------------------------------------------+---------------------------------------------------+-----------------------------+
| :math:`a`                                     | Ice age, ``ice_age`` [days]                       | Thermodynamics, advection   |
+-----------------------------------------------+---------------------------------------------------+-----------------------------+
| :math:`\mathbf{u}_i = (u_i, v_i)`             | Ice velocity, ``Uice`` and ``Vice`` [m/s]         | Dynamics                    |
+-----------------------------------------------+---------------------------------------------------+-----------------------------+
| :math:`\sigma_{xx}, \sigma_{xy}, \sigma_{yy}` | Internal ice stress [N/m]                         | Dynamics                    |
+-----------------------------------------------+---------------------------------------------------+-----------------------------+

The thicknesses are volumes per unit area of the cell, so the ice itself is
:math:`h/A` thick, and the stress is integrated over the thickness of the ice.

Three more fields are diagnosed every step and carried to the next: the
temperature of the ice or snow surface :math:`T_s`, and the temperature and
salinity :math:`T_0`, :math:`S_0` of the thin layer of water in contact with
the underside of the ice. With ``remora.ice.do_momentum = false`` the velocity
and the stress are not solved for and the ice does not move.

The ocean supplies its surface temperature :math:`T_w`, salinity :math:`S_w`,
velocity and free surface, and receives a surface stress, a heat flux and a
salt flux in return.

Thermodynamics
--------------

The thermodynamics treats the ice-covered part of each cell as one layer of
ice of thickness :math:`\hat h_i = h/A + 0.05` m under one layer of snow of
thickness :math:`\hat h_s = h_s/A`, with conductivities :math:`k_i` and
:math:`k_s`. The conductivity of the ice decreases with its brine fraction
:math:`r_b`, which is diagnosed from :math:`T_i` and the ice salinity
:math:`S_i`, the smaller of 3.2 and :math:`S_w`, and it is multiplied by
:math:`2/(1 + e^{-h^2})`, with :math:`h` in meters, which doubles it for thick
ice.

The symbols :math:`\rho_w` and :math:`c_{pw}` are the water density and heat
capacity of the ice module, 1026 kg/m\ :sup:`3` and 3990 J/kg/K unless
``use_rho0`` is set (see :ref:`sec:SeaIce-departures`); :math:`\rho_0` is the
reference density of the run.

**Surface temperature.** The heat the atmosphere delivers to the surface is
linearized in the surface temperature, as :math:`R_a - C_a T_s`, with
:math:`R_a` and :math:`C_a` computed in the bulk flux from the sensible,
latent, shortwave and longwave terms. Setting it equal to the heat conducted
down into the ice gives

.. math::

   T_s = \frac{R_a + b\,T_i}{C_a + b},
   \qquad
   b = \frac{2 k_i}{\hat h_i \left(1 + C_{oa}\right)},
   \qquad
   C_{oa} = \frac{2 k_i \hat h_s}{k_s \hat h_i},

with the temperatures in kelvin and :math:`T_s` kept between -45 C and 0 C.

**Interior temperature.** The ice warms or cools by the difference between
the heat conducted in from the water below and out to the surface above:

.. math::

   \rho_i c_i \hat h_i \frac{\partial T_i}{\partial t}
   = \frac{2 k_i}{\hat h_i}
     \left[ \left(T_0 - T_i\right) + \frac{T_s - T_i}{1 + C_{oa}} \right],
   \qquad
   c_i = c_{pi} + \frac{m\,S_i L_f}{T_i^2},

where :math:`\rho_i` is the ice density, :math:`c_{pi}` the heat capacity of
pure ice, :math:`L_f` the latent heat of fusion, and :math:`m = 0.0543` C per
unit of salinity the slope of the freezing line, :math:`T_f = -m S`. The
second term of :math:`c_i` is the heat taken up by brine pockets.

**Growth and melt.** Ice forms and melts at three interfaces. The rates are
in meters of ice per second; :math:`W_{io}` and :math:`W_{ao}` are positive
for freezing and :math:`W_{ai}` is positive for melting.

.. math::

   W_{io} = \frac{1}{L}\left[\frac{Q_{io}}{\rho_w}
            + c_{pw} C_h \left(T_0 - T_w\right)\right],
   \qquad
   Q_{io} = \frac{2 k_i \left(T_0 - T_i\right)}{\hat h_i},

.. math::

   W_{ao} = \frac{Q_{ao}}{\rho_w L},
   \qquad
   W_{ai} = \frac{b\left(T_i - T_s\right) - Q_{ai}}{\rho_w L}.

:math:`W_{io}` is the growth at the base of the ice: the heat conducted up
into the ice, :math:`Q_{io}`, less the heat the ocean supplies through the
transfer coefficient :math:`C_h`. :math:`W_{ao}` is the growth in open water,
from the heat :math:`Q_{ao}` the open water loses to the atmosphere; it is
applied where the water is at its freezing point and losing heat, or above it
and gaining heat next to ice. :math:`W_{ai}` is the melt at the top of the
ice, from the excess of the atmospheric heat over what the ice conducts away,
where :math:`Q_{ai}` is the heat the ice loses to the atmosphere. The form
given is for bare ice: over a melt pond the density is that of the pond water,
over snow the same heat melts snow instead, and rain that falls on bare ice
freezes onto it. Surface melt removes snow first; snowmelt and melted ice collect
in the melt pond, and what exceeds the smaller of ``max_hmelt`` and :math:`h`
runs off at the rate :math:`W_{ro}`. :math:`L`
is the latent heat of fusion corrected for the brine fraction and for the
heat needed to bring the ice to the temperature of the water (in
:math:`W_{ai}`, of the surface, with :math:`T_s` in place of :math:`T_0`):

.. math::

   L = L_f \left(1 - r_b\right) + c_{pw} T_0
       - \left[\left(1 - r_b\right) c_{pi} + r_b\,c_{pw}\right] T_i .

**Water under the ice.** The layer in contact with the ice is at its freezing
point, :math:`T_0 = -m S_0`, and its salinity balances the salt mixed up from
the ocean against the fresh or salty water the ice releases:

.. math::

   S_0 = \frac{C_s S_w + \left(W_{ai}^{+} - W_{io}\right) S_i}
              {C_s + W_{ai}^{+} + W_{ro} - W_{io}},
   \qquad
   C_{h,s} = \frac{u_\tau}
             {Pr_t \left[\dfrac{1}{\kappa}\ln\dfrac{\Delta z}{z_0} + B_{h,s}\right]},

where :math:`C_h` and :math:`C_s` are the transfer coefficients for heat and
for salt, :math:`W_{ai}^{+} = \max(W_{ai}, 0)`, :math:`u_\tau` is the friction
velocity from the surface stress, :math:`\kappa = 0.4` the von Karman
constant, :math:`Pr_t = 0.85` the turbulent Prandtl number, :math:`\Delta z`
the distance from the surface to the center of the top ocean layer,
:math:`z_0` the roughness of the underside of the ice, and :math:`B_h`,
:math:`B_s` the resistances to heat and to salt of the thin layer next to the
ice in which molecular diffusion dominates.

**Ice volume and concentration.** With the frazil rate :math:`W_{fr}`,

.. math::

   \frac{\partial h}{\partial t}
   = A \left(W_{io} - W_{ai}\right) + \left(1 - A\right) W_{ao} + W_{fr},

.. math::

   \frac{\partial A}{\partial t}
   = \left(1 - A\right)\left(\Phi\,W_{ao} + W_{fr}\right),
   \qquad
   \Phi = \begin{cases} 3 & W_{ao} \ge 0 \\ 0.5 & W_{ao} < 0 \end{cases}

so open water closes faster when it freezes than it opens when it melts. Snow
accumulates on the ice at the snowfall rate times :math:`A`, and the age
advances by the time step, less the share of the ice that is newly formed.

**Frazil ice.** In each ocean layer colder than its freezing point
:math:`T_f = -m S`, the temperature is brought back toward :math:`T_f` and
the heat released becomes ice. The ice formed in the column, less what melts
in warmer layers above, gives :math:`W_{fr}`.

Dynamics
--------

the ice velocity follows a momentum
balance per unit area, with ice mass :math:`M = \rho_i h`:

.. math::

   M \frac{\partial u_i}{\partial t}
   = M f v_i - M g \frac{\partial z_s}{\partial x}
     + \tau_{a}^{x} + A \rho_0 C_{io} \left(u_w - u_i\right)
     + \frac{\partial \sigma_{xx}}{\partial x}
     + \frac{\partial \sigma_{xy}}{\partial y},

.. math::

   M \frac{\partial v_i}{\partial t}
   = - M f u_i - M g \frac{\partial z_s}{\partial y}
     + \tau_{a}^{y} + A \rho_0 C_{io} \left(v_w - v_i\right)
     + \frac{\partial \sigma_{xy}}{\partial x}
     + \frac{\partial \sigma_{yy}}{\partial y}.

The terms on the right are the Coriolis force, the tilt of the sea surface
:math:`z_s`, the wind stress on the ice, the drag of the water, and the
divergence of the internal stress. There is no advection of momentum.

* The wind stress is
  :math:`\boldsymbol{\tau}_a = A \chi \rho_a C_{ai} |\mathbf{U}_a| \mathbf{U}_a`,
  with the drag coefficient ``Cd_ai`` and a factor :math:`\chi` that grows
  with ice thickness (:math:`\chi = 1` with ``mom_bulk``).
* :math:`\mathbf{u}_w` is the ocean velocity averaged over the top 10 m, and
  :math:`C_{io}` is the ice-ocean transfer coefficient, in m/s: by default
  from a logarithmic velocity profile under the ice, or ``Cd_io`` times the
  speed of the ice relative to the water with ``mom_bulk``. In both that speed
  is taken no smaller than 0.15 m/s.
* :math:`f` is the Coriolis parameter and :math:`g` the acceleration of
  gravity.
* In the acceleration and surface-tilt terms :math:`h` is taken no smaller
  than 0.1 m, and in the drag term :math:`A` no smaller than 0.01.

**Rheology.** The internal stress responds to the strain rates

.. math::

   \dot\epsilon_{xx} = \frac{\partial u_i}{\partial x},
   \qquad
   \dot\epsilon_{yy} = \frac{\partial v_i}{\partial y},
   \qquad
   \dot\epsilon_{xy} = \frac{1}{2}\left(\frac{\partial u_i}{\partial y}
                       + \frac{\partial v_i}{\partial x}\right),

through the elastic-viscous-plastic (EVP) law of Hunke and Dukowicz (1997):

.. math::

   \frac{1}{E}\frac{\partial \sigma_{ij}}{\partial t}
   + \frac{\sigma_{ij}}{2\eta}
   + \frac{\eta - \zeta}{4 \eta \zeta}\,\sigma_{kk}\,\delta_{ij}
   + \frac{P}{4\zeta}\,\delta_{ij}
   = \dot\epsilon_{ij} .

Ice behaves as a viscous-plastic material: it resists compression and shear
up to a strength :math:`P`, beyond which it deforms. The elastic term, the
first on the left, is a numerical device that allows the stress to be stepped
explicitly. The strength :math:`P`, the bulk viscosity :math:`\zeta` and the
shear viscosity :math:`\eta` are

.. math::

   P = P^{*} h\, e^{-C\left(1 - A\right)},
   \qquad
   \zeta = \frac{P}{2\Delta},
   \qquad
   \eta = \frac{\zeta}{e^2},

.. math::

   \Delta = \left[\left(\dot\epsilon_{xx} + \dot\epsilon_{yy}\right)^2
            + \frac{1}{e^2}\left(\left(\dot\epsilon_{xx} - \dot\epsilon_{yy}\right)^2
            + 4\,\dot\epsilon_{xy}^2\right)\right]^{1/2},

with ``Pstar`` for :math:`P^{*}`, ``Astrength`` for :math:`C`, ``ellip_sq``
for :math:`e^2`, and :math:`\zeta` kept between ``zetaMin`` and ``zetaMax``.
The elastic parameter is :math:`E = 2 E_0 \rho_i h\,\Delta x^2 / \Delta t_e^2`
with :math:`E_0 = 0.25`. Where the cell-averaged thickness is 0.01 m or less,
the stress is set directly to its viscous-plastic value,

.. math::

   \sigma_{ij} = 2 \eta \dot\epsilon_{ij}
               + \left(\zeta - \eta\right) \dot\epsilon_{kk}\,\delta_{ij}
               - \frac{P}{2}\,\delta_{ij} .

**Subcycling.** The stress and the momentum equations are advanced together
over ``nEVP`` subcycles of length :math:`\Delta t_e` = ``dt/nEVP``. Each
subcycle computes the viscosities from the current velocity, steps the
stress, and steps the velocity with the water drag treated implicitly. The
velocity after the last subcycle is the ice velocity.

The lateral boundary conditions on the three stresses and on the two velocity
components are applied inside every subcycle, to the stress right after it is
stepped and to the velocity right after it is stepped. The previous value a
``radiation`` or ``mixed`` side reads there is the previous subcycle's, so such
a side acts at the elastic step :math:`\Delta t_e`, and the momentum equation
reads the stress the condition has just set outside the domain.

Advection
---------

The ice velocity carries the ice with it. Each of

.. math::

   \phi \in \left\{ A,\; h,\; h_s,\; h_m,\; h\,T_i,\; h\,a \right\}

obeys the conservation law

.. math::

   \frac{\partial \phi}{\partial t} + \nabla \cdot \left(\mathbf{u}_i\,\phi\right) = 0,

solved with MPDATA, an upwind step followed by an antidiffusive correction,
or with the upwind step alone (``remora.ice.advect_upwind``). The temperature and
the age are recovered from the advected products :math:`h T_i` (the enthalpy,
``ice_enthalpy``) and :math:`h a` (``ice_age_thickness``).

What the ocean receives
-----------------------

The ice replaces the surface stress and the surface heat flux of the ocean
and adds to its salt flux. Each is the open-water value weighted by
:math:`1 - A` and the under-ice value weighted by :math:`A`. The ice velocity
in the stress is the one the previous step started from, one step older than
the velocity the momentum equations hold. On the face of a closed side of the
domain the stress is zero.

.. math::

   \frac{\boldsymbol{\tau}_o}{\rho_0}
   = A\,C_{io} \left(\mathbf{u}_i - \mathbf{u}_w\right)
     + \left(1 - A\right) \frac{\boldsymbol{\tau}_{ao}}{\rho_0},

.. math::

   Q_o = \left(1 - A\right) Q_{ao} + A\,Q_{io}
         - \rho_w L \left[A\,W_{io} + \left(1 - A\right) W_{ao}\right],

.. math::

   F_S = \left[A \left(W_{io} - W_{ai}\right) + \left(1 - A\right) W_{ao}
         + W_{fr}\right] \left(S_w - S_i\right)
         - A \left(W_{ro} - W_{ai}^{+}\right) S_w .

:math:`\boldsymbol{\tau}_{ao}` is the wind stress on open water. :math:`Q_o`
is the heat the ocean loses, in W/m\ :sup:`2`: the loss from open water, the
heat conducted into the ice, less the latent heat released by the ice that
forms. :math:`F_S` is the salt flux, in salinity times m/s: freezing leaves
salt behind, because the ice holds a lower salinity :math:`S_i` than the
water, and surface melt and runoff dilute it.

The shortwave radiation :math:`Q_{sw}` that reaches the water is reduced
under the ice-covered part of the cell,

.. math::

   Q_{sw}^{\rm water} = Q_{sw} \left[\left(1 - A\right) + A\,e^{-\tau}\right],

where the optical depth is :math:`\tau = 5\,h/A` for ice up to 0.1 m thick
and :math:`0.5 + 0.1\,(h/A - 0.1)` for thicker ice, plus :math:`20\,h_s/A` for
the snow, with thicknesses in meters.

The options, inputs and outputs of the model are on :ref:`sec:SeaIce`.

References
----------

* Mellor, G. L., and L. Kantha, 1989: An ice-ocean coupled model. *J.
  Geophys. Res.*, **94**, 10937-10954.
* Hunke, E. C., and J. K. Dukowicz, 1997: An elastic-viscous-plastic model
  for sea ice dynamics. *J. Phys. Oceanogr.*, **27**, 1849-1867.
* Ebert, E. E., and J. A. Curry, 1993: An intermediate one-dimensional
  thermodynamic sea ice model for investigating ice-atmosphere interactions. *J.
  Geophys. Res.*, **98**, 10085-10109.
* Budgell, W. P., 2005: Numerical simulation of ice-ocean variability in the
  Barents Sea region. *Ocean Dynamics*, **55**, 370-387.
* Smolarkiewicz, P. K., and W. W. Grabowski, 1990: The multidimensional
  positive definite advection transport algorithm: nonoscillatory option.
  *J. Comput. Phys.*, **86**, 355-375.
