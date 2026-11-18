#include <REMORA.H>

using namespace amrex;

const Vector<REMORA::DVDOutputVar>&
REMORA::dvd_output_vars ()
{
    static const Vector<DVDOutputVar> vars = {
        {"dvd_num_salt", DVD::diag_salt, false,
         "numerical destruction of salinity variance", "second-1"},
        {"dvd_num_temp", DVD::diag_temp, false,
         "numerical destruction of temperature variance", "Celsius2 second-1"},
        {"dvd_phy_salt", DVD::diag_salt, true,
         "destruction of salinity variance by vertical diffusion", "second-1"},
        {"dvd_phy_temp", DVD::diag_temp, true,
         "destruction of temperature variance by vertical diffusion", "Celsius2 second-1"}
    };
    return vars;
}

/**
 * Set the variance decay copies to salinity, temperature and their squares. The copies are
 * then stepped with the tracers, through advection alone, so that at the end of the step they
 * hold A{c} and A{c^2}, with A the advection step (Burchard and Rennau 2008).
 *
 * Both time levels are set, since the predictor reads the previous one, and so is every ghost
 * cell, since the advection stencils read them before the step refills any.
 *
 * Corresponds to the TS_VAR block of pre_step3d.F in the ROMS fork of the beaufort_jet case.
 *
 * @param[in   ] lev    level to operate on
 */
void
REMORA::dvd_reset_copies (int lev)
{
    BL_PROFILE("REMORA::dvd_reset_copies()");
    const int c0 = DVD_comp;
    for (MultiFab* mf : {cons_old[lev], cons_new[lev]}) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box gbx = mfi.growntilebox();
            const Array4<Real> s = mf->array(mfi);
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                const Real salt = s(i,j,k,Salt_comp);
                const Real temp = s(i,j,k,Temp_comp);
                s(i,j,k,c0+DVD::salt)    = salt;
                s(i,j,k,c0+DVD::salt_sq) = salt * salt;
                s(i,j,k,c0+DVD::temp)    = temp;
                s(i,j,k,c0+DVD::temp_sq) = temp * temp;
            });
        }
    }
}

/**
 * Compute the destruction of salinity and temperature variance in the step just taken:
 * the numerical part (A{c^2} - A{c}^2)/dt at cell centers from the copies, and the physical
 * part 2 Akt (dc/dz)^2 at the interior w points from the new tracer and its own vertical
 * diffusivity. Both are zero on land, at the bed and surface w points, and in ghost cells.
 *
 * Corresponds to the TS_VAR block at the end of step3d_t.F in the ROMS fork of the
 * beaufort_jet case, which takes the salinity term with temperature's diffusivity.
 *
 * @param[in   ] lev     level to operate on
 * @param[in   ] dt_lev  time step at this level
 */
void
REMORA::dvd_compute (int lev, Real dt_lev)
{
    BL_PROFILE("REMORA::dvd_compute()");
    const int c0 = DVD_comp;
    const int N = geom[lev].Domain().bigEnd(2);
    const Real odt = one / dt_lev;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_dvd_num[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const Array4<Real const> s    = cons_new[lev]->const_array(mfi);
        const Array4<Real const> mskr = vec_mskr[lev]->const_array(mfi);
        const Array4<Real> num = vec_dvd_num[lev]->array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            const Real m = mskr(i,j,0);
            const Real as = s(i,j,k,c0+DVD::salt);
            const Real at = s(i,j,k,c0+DVD::temp);
            num(i,j,k,DVD::diag_salt) = m * (s(i,j,k,c0+DVD::salt_sq) - as * as) * odt;
            num(i,j,k,DVD::diag_temp) = m * (s(i,j,k,c0+DVD::temp_sq) - at * at) * odt;
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*vec_dvd_phy[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box wbx = mfi.tilebox();
        const Array4<Real const> s    = cons_new[lev]->const_array(mfi);
        const Array4<Real const> Akt  = vec_Akt[lev]->const_array(mfi);
        const Array4<Real const> z_r  = vec_z_r[lev]->const_array(mfi);
        const Array4<Real const> mskr = vec_mskr[lev]->const_array(mfi);
        const Array4<Real> phy = vec_dvd_phy[lev]->array(mfi);

        ParallelFor(wbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // w point k lies between cells k-1 and k; the bed (0) and surface (N+1) have none
            if (k == 0 || k == N+1) {
                phy(i,j,k,DVD::diag_salt) = Real(0.0);
                phy(i,j,k,DVD::diag_temp) = Real(0.0);
                return;
            }
            const Real m = mskr(i,j,0);
            const Real odz = Real(1.0) / (z_r(i,j,k) - z_r(i,j,k-1));
            const Real ds = (s(i,j,k,Salt_comp) - s(i,j,k-1,Salt_comp)) * odz;
            const Real dT = (s(i,j,k,Temp_comp) - s(i,j,k-1,Temp_comp)) * odz;
            phy(i,j,k,DVD::diag_salt) = m * Real(2.0) * Akt(i,j,k,Salt_comp) * ds * ds;
            phy(i,j,k,DVD::diag_temp) = m * Real(2.0) * Akt(i,j,k,Temp_comp) * dT * dT;
        });
    }
}
