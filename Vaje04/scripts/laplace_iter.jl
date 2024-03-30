using Vaje04, Vaje03

using Plots

#Definiramo robni RobniProblemPravokotnik

robni_problem = RobniProblemPravokotnik(
  LaplaceovOperator(2),
  ((0, π), (0, π)),
  [sin, sin, sin, sin])

  #Poiscemo numericno resitev
  U, x, y = resi_iter(robni_problem, 0.1, korak_jacobi)
  #narisemo graf
  surface(x, y, U)

  #Naredimo se animacijo