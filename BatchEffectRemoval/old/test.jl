using Gadfly
using Distributions

#
# X = rand(3,2)
# w = rand(3)
# K = rand(3)
# Binv = rand(2,2)
# b = zeros(2)

n = 50
r = 5
g = 3
m = [50, 50 ,50]
nu2 = 1
sigma2 = 0.1^2

function func_logit(x)
  l = 1./(1+exp(-x))
  l
end


F0 = rand(n,r)
F_list0 = []
C_list0 = []
psi_list0 = []

for j = 1:g
  push!( F_list0,  F0 + rand(n,r) .* sqrt(sigma2));
  push!( C_list0,  rand(r, m[j])* 0.01 .+ 1.0./(1:r).^3    )
end

spy(F0)
spy(F_list0[1])
spy(F_list0[2])
spy(F_list0[3])


function LowerTriange2UpperTriange(x)
  for i=1:size(x,1)
    for j=1:(i-1)
      x[j,i]=x[i,j]
    end
  end
  x
end

#Compute psi = FCF
for j = 1:g
  psi_list_local=[]
  for x = 1: size(C_list0[j],2)
    F_local = F_list0[j]
    push!(psi_list_local, F_local *  ( F_local' .* C_list0[j][:,x]))
  end
  push!(psi_list0, psi_list_local)
end

A_list0= []
for j = 1:g
  A_list_group = []
  for x = 1: size(C_list0[j],2)
    psi_local = psi_list0[j][x]
    A = LowerTriange2UpperTriange(rand(Uniform(0,1),n,n))
    A = (A .< psi_local) .* 1
    push!(A_list_group, A)
  end
  push!(A_list0,A_list_group)
end

spy(A_list0[2][3])


include("randomfactor.jl")

runTest = runBatchRemoval(A_list0, r, 10)

runTest["sigma2"]
plot(y=runTest["trace_sigma2"])


spy(A_list0[1][2])

C_list = runTest["C_list"]
F_list = runTest["F_list"]


plot( x= vec(func_logit(psi_list0[1][3])),
y= vec(func_logit(F_list[1] * diagm(C_list[1][:,3]) * F_list[1]')))


spy(func_logit(F_list[1] * diagm(C_list[1][:,3]) * F_list[1]'))


C_list[1]

spy(F_list[1])
spy(F_list[2])
spy(F_list[3])

spy(runTest["F"])


spy(runTest["C_list"][1])

spy(A_list0[1][3])
spy(func_logit(psi_list0[1][3]))


A_list[1][1]
psi_list[1][1]
