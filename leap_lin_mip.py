import matplotlib.pyplot as plt
import numpy as np
import scipy as sio
import pandas as pd

from dwave.system import LeapHybridCQMSampler
from dimod import ConstrainedQuadraticModel, BinaryQuadraticModel, QuadraticModel

import gurobipy as gp
from gurobipy import GRB
import scipy.io as sio


def populate_and_solve(m):
    data_dir = sio.loadmat('mip.mat')
    f = data_dir['f']

    N = len(f[0])
    rhoTuple = m.addVars(N, vtype=GRB.BINARY)

    objExpr = 0
    for i in range(N):
        objExpr += f[0, i] * rhoTuple[i]

    A = data_dir['obj'][0, 0]
    consExpr = 0
    for i in range(0, N):
        consExpr += rhoTuple[i]
    consExpr -= data_dir['vol'][0, 0]
    consExpr = A * (consExpr)**2

    m.setObjective(objExpr + consExpr, GRB.MINIMIZE)
    m.optimize()

    x = np.zeros((N, 1))
    for i in range(0, N):
        x[i] = rhoTuple[i].X

    return x


data_dir = sio.loadmat('mip.mat')
f = data_dir['f']

A = data_dir['obj'][0, 0]

totalNum = len(f[0])
Q = np.zeros((totalNum, totalNum))
c = np.zeros((totalNum, ))

# objective function
coefficient = np.zeros((totalNum, 1))
for i in range(totalNum):
    coefficient[i] = f[0, i]

c = coefficient

# volume constraint
# coefficient = np.ones((totalNum, 1))
vol = data_dir['vol'][0, 0]
# Q = A * coefficient * coefficient.T
# c += A * (coefficient**2 - vol*coefficient)
# print(Q)
# print(c)

obj = BinaryQuadraticModel(vartype='BINARY')
constraint = QuadraticModel()
for i in range(totalNum):
    obj.add_variable(i)
    constraint.add_variable('BINARY', i)

for i in range(totalNum):
    obj.set_linear(i, c[i])
    constraint.set_linear(i, 1.0)
    # for j in range(i+1, totalNum):
    #     obj.set_quadratic(i, j, Q[i, j] + Q[j, i])

cqm = ConstrainedQuadraticModel()
cqm.set_objective(obj)
cqm.add_constraint(constraint, sense="==", rhs=vol)

sampler = LeapHybridCQMSampler(
    token="DEV-7521f7e94cf42f8fb0430b3d0bfbca3b00264a27")

sampleset = sampler.sample_cqm(cqm, label='lin mip')

feasible_sampleset = sampleset.filter(lambda row: row.is_feasible)

if not len(feasible_sampleset):
    raise ValueError("No feasible solution found")

x = np.zeros((totalNum, 1))
best = feasible_sampleset.first
selected_item_indices = [key for key, val in best.sample.items() if val == 1.0]
x[selected_item_indices] = 1

connection_params = {
    # For Compute Server you need at least this
    #       "ComputeServer": "<server name>",
    #       "UserName": "<user name>",
    #       "ServerPassword": "<password>",

    # For Instant cloud you need at least this
    #       "CloudAccessID": "<access id>",
    #       "CloudSecretKey": "<secret>",
}

with gp.Env(params=connection_params) as env:
    with gp.Model(env=env) as model:
        xCompare = populate_and_solve(model)

print(np.linalg.norm(xCompare - x))
print(best.energy)

print(np.sum(x))

fig = plt.figure()
fig.set_figheight(2)
fig.set_figwidth(6)
plt.imshow(1 - x.reshape(60, 20).T, cmap='gray', vmin=0, vmax=1)
plt.axis("off")
fig.tight_layout()
plt.savefig("leap_lin_mip.eps")
