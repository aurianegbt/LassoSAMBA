H1.all = list(ka=character(),
              V="C1",
              Cl=c("C1","C2"))

param = c("ka","V","Cl")

t.param = c(ka = latex2exp::TeX(r"($ka$)",output="character"),
            V = latex2exp::TeX(r"($V$)",output="character"),
            Cl  = latex2exp::TeX(r"($Cl$)",output="character"))
