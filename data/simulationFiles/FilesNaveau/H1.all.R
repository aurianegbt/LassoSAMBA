H1.all = list(
  phi2 =c("C3","C4","C5"),
  phi1=c("C1","C2","C3"))

param = c("phi2","phi1")

t.param = c(phi2 = latex2exp::TeX(r"($\varphi_2$)",output="character"),
            phi1 = latex2exp::TeX(r"($\varphi_1$)",output="character"))
