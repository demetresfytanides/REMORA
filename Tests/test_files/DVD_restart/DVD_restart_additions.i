# A checkpoint holds no variance decay copies; a restart rebuilds them from the state
remora.do_dvd = true
remora.check_int = 5
remora.plot_int = 10
remora.plot_vars_3d = salt temp x_velocity y_velocity dvd_num_salt dvd_num_temp dvd_phy_salt dvd_phy_temp
