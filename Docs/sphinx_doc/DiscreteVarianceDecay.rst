.. _sec:dvd:

Discrete Variance Decay
=======================

With ``remora.do_dvd = true`` REMORA computes, for salinity and temperature, the destruction of
tracer variance by the discrete advection scheme and by vertical diffusion, following Burchard and
Rennau (2008) and Klingbeil et al. (2014). The diagnostic does not alter the solution.

Variance budget
---------------

Let :math:`c` be salinity or temperature,

.. math::

   \frac{\partial c}{\partial t} + \nabla\cdot(\mathbf{u}\, c) =
   \frac{\partial}{\partial z}\left(\kappa_v \frac{\partial c}{\partial z}\right) + \ldots

In a volume :math:`V` with no boundary fluxes, advection conserves :math:`\int_V c^2\, dV` and
vertical diffusion destroys it:

.. math::

   \frac{d}{dt}\int_V c^2\, dV = -\int_V \chi^{phy}\, dV, \qquad
   \chi^{phy} = 2\,\kappa_v \left(\frac{\partial c}{\partial z}\right)^2 \geq 0.

The discrete advection operator does not conserve :math:`\int_V c^2\, dV`. Over a step
:math:`\Delta t`,

.. math::

   \int_V \left(c^n\right)^2 dV - \int_V \left(c^{n+1}\right)^2 dV =
   \Delta t \int_V \left(\chi^{num} + \chi^{phy}\right) dV,

with :math:`\chi^{num}` the destruction by advection and :math:`\chi^{phy}` that by diffusion.
:math:`\chi^{num}` may be negative.

Numerical destruction
---------------------

Let :math:`\mathcal{A}` denote one advection step, :math:`c^n \mapsto c^{n+1}`, with the model's
scheme, velocities, layer thicknesses and boundary conditions. Then

.. math::

   \chi^{num} = \frac{\mathcal{A}\{c^2\} - \left(\mathcal{A}\{c\}\right)^2}{\Delta t}.

Since :math:`\mathcal{A}` is conservative, :math:`\int_V \mathcal{A}\{c^2\}\, dV = \int_V (c^n)^2\,
dV`, and in the absence of diffusion :math:`\Delta t \int_V \chi^{num}\, dV` equals the loss of
:math:`\int_V c^2\, dV` over the step.

:math:`\mathcal{A}\{c\}` and :math:`\mathcal{A}\{c^2\}` are carried as four internal tracers, copies
of :math:`S`, :math:`S^2`, :math:`T`, :math:`T^2`:

1. At the start of each step the copies are set from the model state, at time levels :math:`n` and
   :math:`n-1` (the predictor reads both) and in all ghost cells.
2. They are advected with the tracers and take the boundary conditions of the tracer they copy. They
   are not diffused and receive no surface or bottom flux, nudging or river concentration.
3. At the end of the tracer step :math:`\chi^{num}` is formed at cell centres.

The reset in step 1 makes :math:`\chi^{num}` a one-step quantity; without it the copies would hold
:math:`\mathcal{A}^k\{c_0^2\}` and :math:`\mathcal{A}^k\{c_0\}` after :math:`k` steps. Diffusion and
forcing are excluded from the copies because they are not part of :math:`\mathcal{A}`; diffusion is
accounted for by :math:`\chi^{phy}`.

Physical destruction
--------------------

At the interior w points between layers :math:`k-1` and :math:`k`,

.. math::

   \chi^{phy}_{k-1/2} = 2\,\kappa_{v,k-1/2}
   \left(\frac{c^{n+1}_k - c^{n+1}_{k-1}}{z_k - z_{k-1}}\right)^2,

with :math:`\kappa_v` the vertical diffusivity of that tracer (``Akt``, background included) and
:math:`z_k` the depth of layer :math:`k`. :math:`\chi^{phy} = 0` at the bed and surface w points.
Horizontal diffusion is not included.

Output
------

+---------------------------+----------------+------------------------------------------------+
| Name                      | Location       | Field                                          |
+===========================+================+================================================+
| ``dvd_num_salt``          | cell centres   | :math:`\chi^{num}` of salinity,                |
|                           |                | second\ :sup:`-1`                              |
+---------------------------+----------------+------------------------------------------------+
| ``dvd_num_temp``          | cell centres   | :math:`\chi^{num}` of temperature,             |
|                           |                | degC\ :sup:`2` second\ :sup:`-1`               |
+---------------------------+----------------+------------------------------------------------+
| ``dvd_phy_salt``          | w points       | :math:`\chi^{phy}` of salinity                 |
+---------------------------+----------------+------------------------------------------------+
| ``dvd_phy_temp``          | w points       | :math:`\chi^{phy}` of temperature              |
+---------------------------+----------------+------------------------------------------------+

The NetCDF history contains all four (``s_rho``, ``s_w``). AMReX plotfiles contain those listed in
``remora.plot_vars_3d``, the w-point fields averaged to cell centres. All four are zero on land and
in the ghost ring written to NetCDF files.

With ``remora.v`` :math:`\geq 1` the integrated-quantity log adds ``salt_var`` and ``temp_var``,
:math:`\int_V (c - \bar c)^2\, dV` with :math:`\bar c` the volume mean, and ``dvd_num_salt``,
``dvd_num_temp``, :math:`\int_V \chi^{num}\, dV`. With advection alone,

.. math::

   \texttt{var}(n-1) - \texttt{var}(n) = \Delta t \int_V \chi^{num}\, dV,

which the regression test ``DVD_budget`` checks.

Restarts and limitations
------------------------

Checkpoints contain :math:`\chi^{num}` and :math:`\chi^{phy}` of the last step but not the copies, which
are rebuilt from the state; restarts are exact, with or without ``remora.do_dvd``.

Inflow and file-driven boundary conditions (``clamped``, ``orlanski_rad_nudg``) for temperature and
salinity are rejected, since no boundary value exists for :math:`c^2`.

Differences from the ``beaufort_jet`` ROMS implementation (``TS_VAR``)
------------------------------------------------------------------------

- :math:`\chi^{phy}` of salinity uses the salinity diffusivity; ROMS uses the temperature diffusivity.
- :math:`\chi^{num}` is stored in separate fields; ROMS stores it in passive-tracer slots, which are
  advected in the following step.
- The copies are reset in all ghost cells; ROMS resets them in the tile interior only.
- The copies are never diffused horizontally; ROMS applies the passive-tracer ``TNU2``.
- Boundary conditions are applied to the copies in tracer units; ROMS applies them to
  thickness-weighted values.
- Only the copies skip vertical diffusion; in ROMS all passive tracers do.

References
----------

- H. Burchard and H. Rennau, 2008: Comparative quantification of physically and numerically induced
  mixing in ocean models. *Ocean Modelling* 20, 293-311.
- K. Klingbeil, M. Mohammadi-Aragh, U. Gräwe and H. Burchard, 2014: Quantification of spurious
  dissipation and mixing -- Discrete variance decay in a Finite-Volume framework. *Ocean Modelling*
  81, 49-64.
