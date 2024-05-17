H1.all = list(V="C1",
              ka=character(),
              Cl=c("C1","C2"))

param = c("V","ka","Cl")

t.param = c(V = latex2exp::TeX(r"($V$)",output="character"),
            ka = latex2exp::TeX(r"($ka$)",output="character"),
            Cl  = latex2exp::TeX(r"($Cl$)",output="character"))
