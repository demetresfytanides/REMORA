.. role:: cpp(code)
  :language: c++

.. _sec:SeaIce:

*******
Sea Ice
*******

REMORA can carry a sea-ice model coupled to the ocean through the surface
stress, heat flux and salt flux. It is a C++ port of the ROMS ice model in
``ROMS/Nonlinear/SeaIce`` (the ``develop`` branch of
`myroms/roms <https://github.com/myroms/roms>`_): the thermodynamics of Mellor
and Kantha (1989), the elastic-viscous-plastic (EVP) dynamics of Hunke and
Dukowicz (1997) as implemented by Budgell (2005), and advection of the ice by
the MPDATA scheme of Smolarkiewicz and Grabowski (1990). Names and defaults
follow the original, so the two can be compared; see `Validation against ROMS`_.
:ref:`sec:SeaIce-equations` gives what the model solves for.

.. warning::

   The sea-ice model is experimental. It runs on a single level only: the ice
   state is not interpolated to, or averaged down from, refined levels, and a
   run with ``amr.max_level > 0`` is refused.

Selecting a model
=================

.. code:: shell

   remora.ice_model = mellor_kantha       # none (the default) disables the ice

+------------------------+----------------------------------------------------------------------+
| Value                  | Model                                                                |
+========================+======================================================================+
| ``none`` (default)     | No ice; identical to a run without the ice model                     |
+------------------------+----------------------------------------------------------------------+
| ``mellor_kantha``      | The ROMS sea-ice model: Mellor-Kantha thermodynamics,                |
|                        | elastic-viscous-plastic (EVP) dynamics and advection                 |
+------------------------+----------------------------------------------------------------------+

Which parts of it run is set by three switches, which correspond to ROMS's
``ICE_THERMO``, ``ICE_MOMENTUM`` and ``ICE_ADVECT`` and are all on by default:
``remora.ice.do_thermo``, ``remora.ice.do_momentum`` and ``remora.ice.do_advect``. Thermodynamic ice
with no velocity is ``remora.ice.do_momentum = false`` with ``remora.ice.do_advect = false``;
see `Runs without advection`_ before using it.

The model is chosen with ``remora.ice_model``, and its settings use the prefix ``remora.ice``.
A minimal set of additions to an ocean inputs file:

.. code:: shell

   remora.bulk_fluxes   = true      # the air-ice fluxes are computed here
   remora.ice_model = mellor_kantha
   remora.plot_vars_2d  = zeta Aice ice_thickness snow_thickness Tice

How a step runs
===============

The ice is advanced once per baroclinic step, after the surface fluxes are
assembled and before the ocean right-hand sides. It reads the ocean at the
start of the step, and the fluxes it writes drive the ocean in the same step.
`Source layout`_ lists the code for each stage.

#. **Bulk flux.** Alongside the open-water fluxes, the bulk flux routine
   computes the air-ice heat flux, the two coefficients the thermodynamics
   uses to solve for the ice surface temperature, the shortwave that reaches
   the water through the ice and snow, the air-ice stress and the snowfall.
#. **Ice-ocean shear.** The speed of the ice relative to the ocean velocity
   averaged over the top 10 m.
#. **Surface stress.** The stress on the ocean becomes the open-water wind
   stress weighted by ``1 - Aice`` plus the ice-water drag weighted by
   ``Aice``, so under full cover the ice-water drag replaces the wind stress.
   This runs whether or not ``remora.ice.do_momentum`` is set.
#. **Dynamics** (``remora.ice.do_momentum``). The elastic-viscous-plastic (EVP) momentum
   equations are advanced over ``nEVP`` elastic subcycles of length
   ``dt/nEVP``, which update the viscosities, the internal stress and the
   ice velocity.
#. **Advection** (``remora.ice.do_advect``). The concentration, the ice, snow and
   melt-pond thicknesses, the enthalpy and the age are advected by the ice
   velocity, then limited to their valid ranges.
#. **Thermodynamics.** Each column grows or melts ice at the air-ice,
   air-ocean and ice-ocean interfaces, updates the snow, the melt ponds and
   the ice temperature, and sets the surface heat and salt fluxes the ocean
   receives. The state is limited again afterwards.
#. **Frazil ice** (optional). At the end of the step, water colder than its
   freezing point is warmed back toward it and the heat is converted to ice. This
   is the only part of the model that changes the three-dimensional tracers.
   It is off by default; see `Departures from ROMS`_.

The ice time step is the baroclinic time step.

State
=====

The prognostic fields are two-dimensional. They use the names ROMS gives them
in ``varinfo.yaml``, which are also the names used for initial files, boundary
files and NetCDF output.

+---------------------------+---------------+------------------------------------------------------+
| Name                      | Location      | Meaning                                              |
+===========================+===============+======================================================+
| ``Aice``                  | rho           | Fraction of the cell covered by ice                  |
+---------------------------+---------------+------------------------------------------------------+
| ``ice_thickness``         | rho           | Ice thickness averaged over the cell [m]             |
+---------------------------+---------------+------------------------------------------------------+
| ``snow_thickness``        | rho           | Snow thickness averaged over the cell [m]            |
+---------------------------+---------------+------------------------------------------------------+
| ``meltpond_thickness``    | rho           | Melt-pond water thickness on the ice [m]             |
+---------------------------+---------------+------------------------------------------------------+
| ``ice_age``               | rho           | Age of the ice [days]                                |
+---------------------------+---------------+------------------------------------------------------+
| ``Tice``                  | rho           | Interior ice temperature [C]                         |
+---------------------------+---------------+------------------------------------------------------+
| ``ice_enthalpy``          | rho           | Enthalpy of the ice and brine [C m]                  |
+---------------------------+---------------+------------------------------------------------------+
| ``ice_age_thickness``     | rho           | Ice thickness times age [m day]                      |
+---------------------------+---------------+------------------------------------------------------+
| ``ice_Sxx``, ``ice_Sxy``, | rho           | Internal ice stress components [N/m]                 |
| ``ice_Syy``               |               |                                                      |
+---------------------------+---------------+------------------------------------------------------+
| ``Uice``, ``Vice``        | u, v          | Ice velocity [m/s]                                   |
+---------------------------+---------------+------------------------------------------------------+

``ice_thickness``, ``snow_thickness`` and ``meltpond_thickness`` are averages
over the whole cell: the thickness of the ice itself is
``ice_thickness / Aice``.

Forcing
=======

The ice needs ``remora.bulk_fluxes = true``, and a run without it is
refused: the air-ice heat flux, the air-ice stress and the snowfall are
computed in the bulk flux routine. The bulk flux inputs are described under
:ref:`list-of-parameters-bulk-fluxes`.

Longwave radiation
  The ice surface takes its longwave in the same form as the ocean: computed
  from air temperature, humidity and cloud fraction by default, from a
  supplied downward longwave with ``remora.longwave_down``, or from a supplied
  net longwave with ``remora.longwave_is_net``.

Shortwave radiation
  The shortwave is attenuated through the ice and snow before it reaches the
  water. With ``remora.shortwave_penetration = true`` it is the attenuated
  shortwave that penetrates the water column.

Albedo
  ``remora.ice.albedo_type`` selects the ice and snow albedo: ``default``
  or ``ec92`` (Ebert and Curry). Where the concentration is at or below
  ``remora.ice.min_ai`` both return ROMS's open-water albedo, the constant 0.06 of
  ``ice_albedo.F``. The shortwave forcing is the net shortwave over open water,
  as for the ocean without ice.

Snowfall
  Precipitation with a negative rate is snow: where ``rain < 0`` its
  magnitude accumulates on the ice as snow of density ``SnowDryRho``, and the
  open water gives up the latent heat of fusion to melt what falls on it.

Evaporation
  Evaporation from the ocean is scaled by the open-water fraction
  ``1 - Aice``.

Initial conditions
==================

The ice starts ice-free unless the initial state says otherwise.

With ``remora.ic_type = netcdf``, REMORA reads the fifteen ice state variables
ROMS reads, and only those the file holds: ``Aice``, ``ice_thickness``,
``meltpond_thickness``, ``snow_thickness``, ``ice_age``, ``ice_Sxx``,
``ice_Sxy``, ``ice_Syy``, ``Tice``, ``ice_enthalpy``, ``ice_age_thickness``,
``Uice``, ``Vice``, ``Uevp`` and ``Vevp``. A variable missing from the file
keeps its ice-free value, and the internal fields are not read by either
model. The last four change no result: the enthalpy and the age thickness are
recomputed from the thickness, temperature and age at the top of every
advection, and ``Uevp`` and ``Vevp`` are reloaded from the ice velocity on the
first EVP subcycle of every step. They are read so that a ROMS initial or
restart file loads unchanged.

With ``remora.ic_type = analytic``, a problem may set the ice state through
``Problem::init_analytic_seaice``; ``SeaIceBox`` does. Problems that do not
start ice-free, and an ocean cold enough to freeze grows its own ice within a
few steps.

Lateral boundary conditions
===========================

Each ice variable takes a condition per side, in the order west, south, east,
north. One entry applies to all four sides.

.. code:: shell

   remora.ice.bc.Aice          = radiation periodic clamped periodic
   remora.ice.bc.ice_thickness = radiation periodic clamped periodic
   remora.ice.bc.Tice          = gradient  periodic clamped periodic
   remora.ice.bc.Uice          = radiation periodic closed  periodic
   remora.ice.bc.Vice          = closed    periodic closed  periodic

``remora.ice.bc_type`` sets one condition for every variable; a
``remora.ice.bc.<name>`` entry overrides it.

The variables are ``Aice``, ``ice_thickness``, ``meltpond_thickness``,
``snow_thickness``, ``ice_age``, ``ice_Sxx``, ``ice_Sxy``, ``ice_Syy``,
``Tice``, ``Uice``, ``Vice``, ``ice_enthalpy`` and ``ice_age_thickness``. The
default is ``closed`` on every side, which suits a closed basin. A condition
on ``ice_enthalpy`` or ``ice_age_thickness`` has no effect: they take their
boundary values from ``Tice`` and ``ice_age``.

Each condition is applied where the field is updated: to an advected field
after its advection, to the thermodynamic fields after the thermodynamics, and
to the stresses and the velocity inside every elastic subcycle (see
:ref:`sec:SeaIce-equations`).

The conditions are the ice model's own, not the ocean's:

``closed``
  For scalars and ``Tice``, zero gradient. For the velocity normal to the
  side, zero flow; for the tangential velocity, ``bc_gamma2`` times the
  interior value.

``gradient``
  Zero gradient. For ``Tice`` it is applied to the enthalpy, and the
  temperature is recovered from it.

``clamped``
  The boundary value. For ``Tice`` the enthalpy is clamped to the boundary
  ``ice_thickness`` times the boundary ``Tice`` and the temperature is recovered
  from it, so it needs boundary values for both and ``ice_thickness`` clamped
  on that side.

``mixed``
  The boundary value on inflow, an upstream value on outflow. For scalars and
  ``Tice``; on ``Uice`` or ``Vice`` it is rejected, as it applies nothing there.

``radiation``
  A radiation condition, which lets disturbances leave through the side,
  with nudging to the boundary value on the sides named in ``bc_nudging``.
  On ``Tice`` it applies nothing, and REMORA prints a warning.

``periodic``
  A periodic side; nothing is applied.

Boundary values come from a NetCDF file or from constants:

.. code:: shell

   remora.ice.nc_bdry_file            = ice_bry.nc
   remora.ice.bc_value.Aice           = 0.0 0.0 0.5 0.0   # west south east north
   remora.ice.bc_value.ice_thickness  = 0.0 0.0 0.2 0.0

A ``clamped`` or ``mixed`` side, or a ``radiation`` side with nudging, needs a
value. With ``remora.ice.nc_bdry_file`` set, REMORA reads
``<name>_west``, ``<name>_east``, ``<name>_south`` and ``<name>_north`` for
the sides that need them, and interpolates in time. A variable whose records
are not all in the file takes ``remora.ice.bc_value.<name>`` instead,
which defaults to zero; REMORA prints which source each variable uses. The
time coordinate is ``Aice_time`` for ``Aice`` and ``ocean_time`` for every
other variable, and can be changed per variable with
``remora.ice.bc_time_name.<name>``. Reading a file requires a PnetCDF
build, and the file needs at least two time records.

+--------------------------------------+-------------+---------------------------------------------+
| Option                               | Default     | Meaning                                     |
+======================================+=============+=============================================+
| ``remora.ice.bc_nudging``            | 0           | Nudge a ``radiation`` side; 1 or 4 entries  |
+--------------------------------------+-------------+---------------------------------------------+
| ``remora.ice.bc_tau_in``             | 0.0         | Nudging rate on inflow [1/s]                |
+--------------------------------------+-------------+---------------------------------------------+
| ``remora.ice.bc_tau_out``            | 0.0         | Nudging rate on outflow [1/s]               |
+--------------------------------------+-------------+---------------------------------------------+
| ``remora.ice.bc_gamma2``             | 1.0         | Tangential slip: 1 free slip, -1 no slip    |
+--------------------------------------+-------------+---------------------------------------------+

``bc_gamma2`` is the slip of the tangential ice velocity at a ``closed`` side
and is separate from the ocean's slip condition.

Output
======

AMReX plotfiles
---------------

Ice fields are requested by name in ``remora.plot_vars_2d``:

.. code:: shell

   remora.plot_vars_2d = zeta Aice ice_thickness snow_thickness Tice ice_sst

The names available are ``Aice``, ``ice_thickness``, ``snow_thickness``,
``meltpond_thickness``, ``ice_age``, ``Tice``,
``ice_sst`` (the temperature of the ice or snow surface), ``under_ice_temp``
and ``under_ice_salt`` (the temperature and salinity of the molecular
sublayer under the ice) and ``io_massflux`` (the ice-ocean mass flux).
``seaice`` stands for all of them. They are available only when the ice model is
on.

NetCDF history
--------------

With ``remora.plotfile_type = netcdf`` and the ice model on, the history file
carries the whole ice state without being asked. It holds the fields of
`State`_ and:

+--------------------------------------+--------------------------------------------------------+
| Name                                 | Meaning                                                |
+======================================+========================================================+
| ``ice_sst``                          | Temperature of the ice or snow surface [C]             |
+--------------------------------------+--------------------------------------------------------+
| ``under_ice_temp``                   | Temperature of the sublayer under the ice [C]          |
+--------------------------------------+--------------------------------------------------------+
| ``under_ice_salt``                   | Salinity of the sublayer under the ice                 |
+--------------------------------------+--------------------------------------------------------+
| ``io_friction_vel``                  | Ice-ocean friction velocity [m/s]                      |
+--------------------------------------+--------------------------------------------------------+
| ``io_massflux``                      | Ice-ocean mass flux [m/s]                              |
+--------------------------------------+--------------------------------------------------------+
| ``io_transfer_coef``                 | Ice-ocean momentum transfer coefficient [m/s]          |
+--------------------------------------+--------------------------------------------------------+
| ``ice_snow_heat_conductivity``       | Conductivity of the ice and snow [W/m2/K]              |
+--------------------------------------+--------------------------------------------------------+
| ``ice_snow_heat_flux``               | Heat-balance right-hand side over ice [W/m2]           |
+--------------------------------------+--------------------------------------------------------+
| ``ice_divergence``                   | Thickness change by advection per unit time [m/s]      |
+--------------------------------------+--------------------------------------------------------+
| ``ai_melt_freeze``                   | Melt or freeze rate, air-ice interface [m/s]           |
+--------------------------------------+--------------------------------------------------------+
| ``ao_melt_freeze``                   | Melt or freeze rate, air-ocean interface [m/s]         |
+--------------------------------------+--------------------------------------------------------+
| ``io_melt_freeze``                   | Melt or freeze rate, ice-ocean interface [m/s]         |
+--------------------------------------+--------------------------------------------------------+
| ``frazil_freeze``                    | Rate of ice accretion by frazil [m/s]                  |
+--------------------------------------+--------------------------------------------------------+
| ``Uevp``, ``Vevp``                   | Ice velocity after the last elastic subcycle [m/s]     |
+--------------------------------------+--------------------------------------------------------+
| ``shflux``                           | Net heat flux into the ocean [W/m2]                    |
+--------------------------------------+--------------------------------------------------------+
| ``stflx_salt``                       | Salt flux into the ocean [salinity m/s]                |
+--------------------------------------+--------------------------------------------------------+
| ``swrad_water``                      | Shortwave reaching the water [W/m2]                    |
+--------------------------------------+--------------------------------------------------------+

``shflux`` and ``stflx_salt`` are the fluxes the ocean received after the ice
replaced them; ``shflux`` includes the shortwave. (In an AMReX plotfile ``shflux``
is the sensible heat flux, and the ``output_forcing`` fields hold the fluxes before
the ice replaces them.) ``ice_divergence`` is
positive where the ice converges, and ``Uevp`` and ``Vevp`` equal ``Uice`` and
``Vice``.

NetCDF output requires a PnetCDF build (see :ref:`sec:Plotfiles`).

Restarts
========

Checkpoints hold the full ice state, both time levels and every internal
field, and a restart continues the ice exactly. The ``Header`` of a checkpoint
written with ice records the ice model:

* A checkpoint written with ice restarts with ice when ``remora.ice_model`` is
  set. Otherwise its ice fields are ignored.
* A checkpoint written without ice, including one from before the ice model
  existed, can be restarted with the ice on. The ocean continues and the ice
  starts ice-free, which REMORA reports. This is the way to switch the ice on
  after an ocean spin-up.

Options
=======

These select among formulations that ROMS chooses at compile time; the effect
names the C-preprocessor option each corresponds to.

+--------------------------------+-------------+----------------------------------------------------+
| Option                         | Default     | Effect                                             |
+================================+=============+====================================================+
| ``remora.ice.do_thermo``       | ``true``    | Run the Mellor-Kantha thermodynamics. With it off  |
|                                |             | the ocean keeps the surface fluxes the bulk flux   |
|                                |             | computed. ROMS: ``ICE_THERMO``.                    |
+--------------------------------+-------------+----------------------------------------------------+
| ``remora.ice.do_momentum``     | ``true``    | Run the EVP rheology and the ice momentum          |
|                                |             | equations. With it off the ice velocity stays as   |
|                                |             | it was initialised. ROMS: ``ICE_MOMENTUM``.        |
+--------------------------------+-------------+----------------------------------------------------+
| ``remora.ice.do_advect``       | ``true``    | Advect the ice. ROMS: ``ICE_ADVECT``.              |
+--------------------------------+-------------+----------------------------------------------------+
| ``remora.ice.advect_upwind``   | ``false``   | Donor-cell upwind advection, without the           |
|                                |             | antidiffusive pass of MPDATA. ROMS:                |
|                                |             | ``ICE_UPWIND``.                                    |
+--------------------------------+-------------+----------------------------------------------------+
| ``remora.ice.mom_bulk``        | ``false``   | Bulk ice-ocean drag: the transfer coefficient is   |
|                                |             | ``Cd_io`` times the ice-ocean shear, in place of   |
|                                |             | the log-layer form that depends on ice thickness.  |
|                                |             | ROMS: ``ICE_MOM_BULK``.                            |
+--------------------------------+-------------+----------------------------------------------------+
| ``remora.ice.convsnow``        | ``false``   | Convert to ice the snow whose base its own weight  |
|                                |             | has pushed below the water line. ROMS:             |
|                                |             | ``ICE_CONVSNOW``.                                  |
+--------------------------------+-------------+----------------------------------------------------+
| ``remora.ice.albedo_type``     | ``default`` | ``ec92`` selects the Ebert and Curry albedo. ROMS: |
|                                |             | ``ICE_ALB_EC92``.                                  |
+--------------------------------+-------------+----------------------------------------------------+
| ``remora.ice.strength_quad``   | ``false``   | Ice strength from a gridded parameter and the      |
|                                |             | square of the thickness. The parameter is never    |
|                                |             | filled, so the ice has no strength; the option     |
|                                |             | exists so that a run built with it can be matched. |
|                                |             | ROMS: ``ICE_STRENGTH_QUAD``.                       |
+--------------------------------+-------------+----------------------------------------------------+

Not ported: ``ICE_SHOREFAST``, ``ICE_KEEL``, ``ICE_BIO``, ``ICESHELF``,
``ALBEDO_CSIM``, ``FASTICE_CLIMATOLOGY``, ``OUTFLOW_MASK``, ``RADIATION_2D``, the
``WET_DRY`` factors, the climatology nudging of the ice, and the mixed-layer depth of
KPP, so that the ice-ocean shear averages over 10 m.

Parameters
==========

Each parameter is set as ``remora.ice.<name>``. Names and defaults are
those of ``ROMS/External/ice.in``.

+----------------------+---------------+------------------------------------------------------------+
| Parameter            | Default       | Meaning                                                    |
+======================+===============+============================================================+
| ``nEVP``             | 60            | Elastic subcycles per ice time step                        |
+----------------------+---------------+------------------------------------------------------------+
| ``IceRho``           | 900.0         | Ice density [kg/m3]                                        |
+----------------------+---------------+------------------------------------------------------------+
| ``SnowDryRho``       | 330.0         | Dry snow density [kg/m3]                                   |
+----------------------+---------------+------------------------------------------------------------+
| ``SnowWetRho``       | 450.0         | Wet snow density [kg/m3]                                   |
+----------------------+---------------+------------------------------------------------------------+
| ``Cd_ai``            | 3.0e-3        | Air-ice drag coefficient                                   |
+----------------------+---------------+------------------------------------------------------------+
| ``Cd_io``            | 7.5e-3        | Ice-ocean drag coefficient                                 |
+----------------------+---------------+------------------------------------------------------------+
| ``Pstar``            | 25000.0       | Ice strength scale [N/m2]                                  |
+----------------------+---------------+------------------------------------------------------------+
| ``Astrength``        | 20.0          | Exponential weighting of the ice strength on concentration |
+----------------------+---------------+------------------------------------------------------------+
| ``zetaMin``          | 0.0           | Lower bound of the bulk viscosity                          |
+----------------------+---------------+------------------------------------------------------------+
| ``zetaMax``          | 1.0e14        | Upper bound of the bulk viscosity                          |
+----------------------+---------------+------------------------------------------------------------+
| ``ellip_sq``         | 4.0           | Squared ellipticity of the yield curve                     |
+----------------------+---------------+------------------------------------------------------------+
| ``min_ai``           | 0.0           | Minimum ice concentration                                  |
+----------------------+---------------+------------------------------------------------------------+
| ``max_ai``           | 1.0           | Maximum ice concentration                                  |
+----------------------+---------------+------------------------------------------------------------+
| ``min_hi``           | 0.0           | Minimum ice thickness [m]                                  |
+----------------------+---------------+------------------------------------------------------------+
| ``max_hmelt``        | 0.10          | Maximum melt-pond thickness [m]                            |
+----------------------+---------------+------------------------------------------------------------+

``min_ai`` must not exceed ``max_ai``, and ``nEVP`` must be at least 1. The
physical constants of the thermodynamics (conductivities, heat capacities,
latent heat of fusion, the freezing-line slope ``-0.0543`` C per unit
salinity) are fixed; they are listed in ``Source/SeaIce/REMORA_SeaIce.H``.

Six more parameters of ``ice.in`` are read and have no effect, in REMORA as in
ROMS: ``AirRho``, ``stressAng``, ``ice_emiss``, ``spec_heat_air``,
``trans_coeff`` and ``sublimation``. The air-ice fluxes take their constants
from the bulk flux (a specific heat of air of 1004.67 J/kg/K, an emissivity of
0.97, a transfer coefficient of 1.75e-3 and a latent heat of sublimation of
2.834e6 J/kg), and no turning angle is applied to the ice-water drag. REMORA
prints a warning when an inputs file sets one of them.

.. _sec:SeaIce-departures:

Departures from ROMS
====================

The port copies the original, including behavior that looks unintended. Where
a more consistent alternative exists, it is one of the options below. Each
defaults to ``false``, which keeps the original behavior.

+------------------------------------------+-----------+----------------------------------------------+
| Option                                   | Default   | Effect of ``true``                           |
+==========================================+===========+==============================================+
| ``remora.ice.do_frazil``                 | ``false`` | Runs the frazil-ice calculation, which ROMS  |
|                                          |           | contains but never calls. Without it ice     |
|                                          |           | forms only through the heat the surface      |
|                                          |           | loses.                                       |
+------------------------------------------+-----------+----------------------------------------------+
| ``remora.ice.pressure_freezing``         | ``false`` | Includes the pressure term in the frazil     |
|                                          |           | freezing point, in place of ``-0.0543 S`` at |
|                                          |           | every depth.                                 |
+------------------------------------------+-----------+----------------------------------------------+
| ``remora.ice.use_rho0``                  | ``false`` | Uses the run's ``remora.rho0`` and heat      |
|                                          |           | capacity in the ice module, in place of a    |
|                                          |           | fixed 1026 kg/m3 and 3990 J/kg/K. In fresh   |
|                                          |           | water with ``rho0 = 1000`` this changes      |
|                                          |           | every ice-derived heat flux by 2.7%.         |
+------------------------------------------+-----------+----------------------------------------------+
| ``remora.ice.albedo_melt_fix``           | ``false`` | Compares the surface temperature with the    |
|                                          |           | melting point in Celsius in the default      |
|                                          |           | albedo. Otherwise a Celsius temperature is   |
|                                          |           | compared with a Kelvin threshold and the     |
|                                          |           | melting branch is never taken.               |
+------------------------------------------+-----------+----------------------------------------------+
| ``remora.ice.zmax_clamp``                | ``false`` | Limits the bulk viscosity to ``2.5e8`` times |
|                                          |           | the ice strength, in addition to             |
|                                          |           | ``zetaMax``.                                 |
+------------------------------------------+-----------+----------------------------------------------+
| ``remora.ice.qnet_ao_open_water_sw``     | ``false`` | Forms the open-water heat flux from the      |
|                                          |           | shortwave that open water receives, in place |
|                                          |           | of the shortwave already attenuated by the   |
|                                          |           | ice cover.                                   |
+------------------------------------------+-----------+----------------------------------------------+
| ``remora.ice.bc_fix.tibc_mixed``         | ``false`` | Applies the ``mixed`` condition on ``Tice``  |
|                                          |           | to an east or south inflow as the west and   |
|                                          |           | north sides apply it.                        |
+------------------------------------------+-----------+----------------------------------------------+
| ``remora.ice.bc_fix.uibc_north_nudging`` | ``false`` | Nudges ``Uice`` on the north side, as on the |
|                                          |           | other three.                                 |
+------------------------------------------+-----------+----------------------------------------------+
| ``remora.ice.bc_fix.vibc_west_closed``   | ``false`` | Scales the interior value in the west        |
|                                          |           | ``closed`` condition on ``Vice``, as the     |
|                                          |           | other three sides do, so ``bc_gamma2`` means |
|                                          |           | the same there.                              |
+------------------------------------------+-----------+----------------------------------------------+
| ``remora.ice.bc_fix.vibc_north_mask``    | ``false`` | Masks the north ``clamped`` and ``gradient`` |
|                                          |           | conditions on ``Vice`` at the boundary       |
|                                          |           | point.                                       |
+------------------------------------------+-----------+----------------------------------------------+

``remora.ice.bc_fix_all = true`` sets the four ``bc_fix`` options
together; one set individually still takes precedence.

Always different from the original:

* The southern ``radiation`` condition on ``Vice`` follows its own entry. ROMS tests the
  condition of another variable there.
* The boundary values are interpolated without rounding the weights to the millisecond.
* The history file holds the internal stress of the last elastic subcycle
  (see `Validation against ROMS`_).
* The ghost cells of the surface temperature, ``Qcon``, ``Qrhs`` and the surface fluxes on an
  open or closed side are not filled. Nothing reads them.

.. _sec:SeaIce-no-advection:

Runs without advection
----------------------

``remora.ice.no_advect_two_states`` (default ``true``) applies when the
ice is not advected, that is with ``remora.ice.do_advect = false``.

ROMS keeps two copies of the ice state and switches which one is current at
the start of every ice step. The advection is what fills the new current copy
from the old one. Built without ``ICE_ADVECT``, it therefore advances the
two copies on alternate steps, each from where it was two steps earlier, and
computes the surface fluxes of one step from the other copy. An initial state
read from a file goes into one copy only. Ice then grows at about half the
rate of a single state, and the surface heat flux alternates from step to
step.

The default reproduces this so that such a run can be matched, and REMORA
prints a warning at startup.

.. warning::

   For a thermodynamics-only run that is not being compared with ROMS, set
   ``remora.ice.no_advect_two_states = false``. The ice then has one
   state.

Validation against ROMS
=======================

To reproduce a ROMS ice run, transcribe ``ice.in`` into the parameters above,
set the options its application defines, and set these two, which REMORA
defaults differently.

.. code:: shell

   remora.g                = 9.81      # REMORA defaults to 9.80665
   remora.bulk_flux_scheme = coare35   # unless the application defines COARE_30

The ``LAKE_ICE`` application defines ``ICE_MOMENTUM``, ``ICE_EVP``,
``ICE_ADVECT``, ``ICE_UPWIND``, ``ICE_MOM_BULK``, ``ICE_CONVSNOW``,
``ICE_ALB_EC92`` and ``SOLAR_SOURCE``, and computes its longwave (``LONGWAVE``), which corresponds
to:

.. code:: shell

   remora.ice_model              = mellor_kantha
   remora.ice.mom_bulk           = true
   remora.ice.convsnow           = true
   remora.ice.advect_upwind      = true
   remora.ice.albedo_type        = ec92
   remora.shortwave_penetration  = true
   remora.jerlov_type            = I

Three things to know when comparing history files:

* A ROMS history record holds the ocean at its own time but the ice and the
  surface fluxes of the following step. Its record ``k`` compares with REMORA
  record ``k + 1`` for the ice and the fluxes, and with REMORA record ``k``
  for the ocean.
* With an even ``nEVP``, on odd steps the ROMS file holds the internal stress (``ice_Sxx``,
  ``ice_Sxy``, ``ice_Syy``) of the next-to-last elastic subcycle. REMORA
  writes the last, which is the one both models carry forward.
* The ocean's net longwave differs by about 7e-5 in relative terms, because
  REMORA's ocean bulk flux uses a more precise Stefan-Boltzmann constant. The
  ice block uses the original value.

After one step of a ``LAKE_ICE`` run started from the same state, the ice
state, the ice velocity and the surface stress and salt flux agree to the
single precision of the ROMS history file, about 6e-8 in relative terms.
Later steps differ as far as the two ocean models do.

Diagnostics
===========

``remora.ice.probe_i`` and ``remora.ice.probe_j`` name one cell whose
ice column is printed every ``remora.ice.probe_int`` steps: the inputs the
thermodynamics received (``ICEPROBE ... IN``), the state it produced
(``ICEPROBE ... OUT``) and the ice velocity (``ICEVEL``).
``remora.ice.probe_bc = true`` also prints the ghost values the lateral
boundary conditions wrote on the four sides at that row and column
(``ICEBC``). The probe is off by default.

Regression coverage
===================

The sea-ice regression tests (27, or 29 in PnetCDF builds), what each asserts, and what no test covers are listed under
:ref:`sec:ci-seaice-tests`.

Source layout
=============

The model is in ``Source/SeaIce``. Kernels that act on one cell or one column
are inline functions in the headers and are called from the drivers.

+-----------------------------+------------------------------------+-------------------------------+
| Contents                    | REMORA files                       | ROMS source                   |
+=============================+====================================+===============================+
| State, parameters, names    | | ``REMORA_SeaIce.H``              | | ``ice_mod.h``               |
|                             | | ``REMORA_SeaIce.cpp``            | | ``ice.in``                  |
+-----------------------------+------------------------------------+-------------------------------+
| Driver for one ice step     | ``REMORA_advance_sea_ice.cpp``     | ``seaice.F``                  |
+-----------------------------+------------------------------------+-------------------------------+
| Air-ice heat flux, albedo   | ``REMORA_IceBulkFlux.H``           | | ``bulk_flux.F``             |
|                             |                                    | | ``ice_albedo.F``            |
+-----------------------------+------------------------------------+-------------------------------+
| Ice-ocean shear and stress  | | ``REMORA_IceStress.H``           | | ``ice_spdiw.F``             |
|                             | | ``REMORA_ice_stress.cpp``        | | ``ice_vbc.F``               |
+-----------------------------+------------------------------------+-------------------------------+
| EVP dynamics                | | ``REMORA_IceEVP.H``              | | ``ice_evp.F``               |
|                             | | ``REMORA_IceElastic.H``          | | ``ice_evp_sig.F``           |
|                             | | ``REMORA_ice_dynamics.cpp``      | | ``ice_elastic.F``           |
+-----------------------------+------------------------------------+-------------------------------+
| Advection, enthalpy,        | | ``REMORA_IceAdvect.H``           | | ``ice_advect.F``            |
| limiter                     | | ``REMORA_ice_advect.cpp``        | | ``ice_smolar.h``            |
|                             |                                    | | ``ice_enthalpy.F``          |
|                             |                                    | | ``ice_limit.F``             |
+-----------------------------+------------------------------------+-------------------------------+
| Thermodynamics of a column  | ``REMORA_IceThermo.H``             | ``ice_mk.h``                  |
+-----------------------------+------------------------------------+-------------------------------+
| Frazil ice                  | | ``REMORA_IceFrazil.H``           | ``ice_frazil.F``              |
|                             | | ``REMORA_ice_frazil.cpp``        |                               |
+-----------------------------+------------------------------------+-------------------------------+
| Lateral boundary conditions | | ``REMORA_IceBC.H``               | | ``ice_bc2d.F``              |
|                             | | ``REMORA_ice_bc.cpp``            | | ``ice_tibc.F``              |
|                             |                                    | | ``ice_uibc.F``              |
|                             |                                    | | ``ice_vibc.F``              |
+-----------------------------+------------------------------------+-------------------------------+
| Boundary values from NetCDF | ``REMORA_ice_bdry_data.cpp``       | ``set_data.F``                |
+-----------------------------+------------------------------------+-------------------------------+
| Initial state from NetCDF   | ``REMORA_ice_init.cpp``            | ``get_state.F``               |
+-----------------------------+------------------------------------+-------------------------------+
