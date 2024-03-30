module Vaje04



function korak_jacobi(L::LaplaceovOperator{2}, U0)
  U = copy(U0)
  n, m = size(U)
  #zanka po notranjih elementih U
  for i=2:n-1, j=2:m-1
    U[i,j] = (U0[i-1, j] + U0[i, j-1] + U0[i+1, j] + U0[i, j+1])/4
  end
  return U
end


function korak_gs(L::LaplaceovOperator{2}, U0)
  U = copy(U0)
  n, m = size(U)
  #zanka po notranjih elementih U
  for i=2:n-1, j=2:m-1
    U[i,j] = (U[i-1, j] + U[i, j-1] + U[i+1, j] + U[i, j+1])/4
  end
  return U
end



function korak_sor(L::LaplaceovOperator{2}, U0, ω)
  U = copy(U0)
  n, m = size(U)
  #zanka po notranjih elementih U
  for i=2:n-1, j=2:m-1
    U[i,j] = (1-ω)*U0[i,j] + ω*(U[i-1, j] + U[i, j-1] + U[i+1, j] + U[i, j+1])/4
  end
  return U
end


function iteracija(korak, U0; maxit=1000, tol=1e-10)
  for i=1:maxit
    U = korak(U0)
    if norm(U - U0) < tol
      #zaporedje je skonvergiralo
      return U, i
    end
    U0 = U
  end
  throw("Zapoeredje ni konvergentno")
end

function resi_iter(robni_problem::RobniProblemPravokotnik, h, korak)
    # poisci zacetni priblizek
    U0, x, y = zacetni_priblizek(robni_problem, h)
    # pozeni iteracijo
    k(U) = korak(robni_problem.operator, U)
    U, it = iteracija(k, U0)
    # vrni rezultat
    U, x, y
end

function zacetni_priblizek(robni_problem::RobniProblemPravokotnik, h)
  ((a,b),(c,d)) = robni_problem.meje
  f_s, f_d, f_z, f_l = robni_problem.rp
  n = Int(floor((b - a) / h)) - 1 #stevilo notranjih vozlisc v x smeri
  m = Int(floor((d - c) / h)) - 1 #stevilo notranjih vozlisc v y smeri
  U = zeros(n + 2, m + 2)
  x = LinRange(a, b, n + 2)
  y = LinRange(c, d, m + 2)
  U[1, :] = f_l.(y)
  U[n + 2, :] = f_d.(y)
  U[:, 1] = f_s.(x)
  U[:, m + 2] = f_z.(x)
  return U, x, y
end



end # module Vaje04
