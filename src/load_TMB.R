TMB::compile("src/TMB/tau2_rho_f_v2.cpp",CPPFLAGS="-Wno-ignored-attributes")
TMB::compile("src/TMB/tau2_rho_f_v1.cpp",CPPFLAGS="-Wno-ignored-attributes")
dyn.load(dynlib("tau2_rho_f_v2"))
dyn.load(dynlib("tau2_rho_f_v1"))