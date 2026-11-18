#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi((mf_Akv), TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
      Array4<Real> const& Akv = (mf_Akv).array(mfi);
      Array4<Real> const& Akt = (mf_Akt).array(mfi);
      Box bx = mfi.growntilebox(IntVect(NGROW,NGROW,0));
      Gpu::streamSynchronize();
      // Constant at the background values, as ROMS's double gyre takes AKV_BAK and AKT_BAK
      const Real Akv_bak = m_solverChoice.Akv_bak;
      const Real Akt_bak_temp = m_solverChoice.Akt_bak[Temp_comp];
      const Real Akt_bak_salt = m_solverChoice.Akt_bak[Salt_comp];
      ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
      {
        Akv(i,j,k) = Akv_bak;

        Akt(i,j,k,Temp_comp) = Akt_bak_temp;
        Akt(i,j,k,Salt_comp) = Akt_bak_salt;
      });
    }
