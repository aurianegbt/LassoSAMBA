H1.all = list(phi_S=c("G1","G2"),
              phi_L=c("G3","G4","G5","G6","G7","G8","G9","G10","G11","G13","G14","G15","G16"),
              delta_AB="G15")

param = c("phi_S","phi_L","delta_AB")

t.param = c(phi_S = latex2exp::TeX(r"($\varphi_S$)",output="character"),
            phi_L = latex2exp::TeX(r"($\varphi_L$)",output="character"),
            delta_AB  = latex2exp::TeX(r"($\delta_{Ab}$)",output="character"))
