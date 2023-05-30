import matplotlib.pyplot as plt
import numpy as np
import scipy as sio
import pandas as pd

from dwave.system import DWaveSampler, EmbeddingComposite
from dwave.system import LeapHybridCQMSampler
from dwave.samplers import SimulatedAnnealingSampler
from dimod import ConstrainedQuadraticModel, BinaryQuadraticModel, QuadraticModel

import gurobipy as gp
from gurobipy import GRB
import scipy.io as sio

from dwave.cloud.sw import Client


def populate_and_solve(m):
    data_dir = sio.loadmat('cons_mip.mat')
    weight = np.array(data_dir['targetWeight'])
    b = np.array(data_dir['targetB'])
    Upper = np.max(np.abs(b))

    N = weight.shape[1]
    rhoTuple = m.addVars(N, vtype=GRB.BINARY)

    M = 10
    etaTuple = m.addVars(M+1, vtype=GRB.BINARY)
    eta = Upper*(1 - 1 / (2**M))*etaTuple[0]
    for i in range(1, M+1):
        eta += Upper/(2**i)*etaTuple[i]

    Nc = weight.shape[0]
    alphaTuple = m.addVars(Nc, M+1, vtype=GRB.BINARY)

    objExpr = eta

    A = Upper
    consExpr = 0
    for c in range(Nc):
        subConsExpr = 0
        for i in range(N):
            subConsExpr += -weight[c, i] * rhoTuple[i]
        alpha = Upper*(1 - 1 / (2**M))*alphaTuple[c, 0]
        for i in range(1, M+1):
            alpha += Upper/(2**i)*alphaTuple[c, i]
        subConsExpr = subConsExpr - eta + alpha - b[c]
        consExpr += A * (subConsExpr)**2

    subConsExpr = 0
    for i in range(0, N):
        subConsExpr += rhoTuple[i] / N
    subConsExpr -= data_dir['targetVol'][0, 0] / N
    consExpr += A * (subConsExpr)**2

    m.setObjective(objExpr + consExpr, GRB.MINIMIZE)
    m.optimize()

    eta = Upper*(1 - 1 / (2**N))*etaTuple[0].X
    for i in range(1, M+1):
        eta += Upper/(2**i)*etaTuple[i].X

    x = np.zeros((N, 1))
    for i in range(0, N):
        x[i] = rhoTuple[i].X

    return x


data_dir = sio.loadmat('cons_mip.mat')
weight = np.array(data_dir['targetWeight'])
b = np.array(data_dir['targetB'])
Upper = np.max(np.abs(b))

N = weight.shape[1]
Nc = weight.shape[0]
M = 10
A = Upper

totalNum = N + M+1 + Nc * (M+1)
Q = np.zeros((totalNum, totalNum))
c = np.zeros((totalNum, ))

etaCoefficient = np.zeros((M+1, ))
etaCoefficient[0] = Upper*(1 - 1 / (2**M))
for i in range(1, M+1):
    etaCoefficient[i] = Upper/(2**i)

alphaCoefficient = np.zeros((M+1, Nc))
for i in range(Nc):
    alphaCoefficient[0, i] = Upper*(1 - 1 / (2**M))
    for j in range(1, M+1):
        alphaCoefficient[j, i] = Upper/(2**j)

obj = BinaryQuadraticModel(vartype='BINARY')
for i in range(totalNum):
    obj.add_variable(i)

for i in range(totalNum):
    obj.set_linear(i, c[i])

# objective function
coefficient = np.zeros((totalNum, 1))
for i in range(M+1):
    coefficient[N+i] = etaCoefficient[i]

for i in range(totalNum):
    obj.set_linear(i, coefficient[i])

cqm = ConstrainedQuadraticModel()
cqm.set_objective(obj)

# cuts
for j in range(Nc):
    coefficient = np.zeros((totalNum, 1))
    constraint = BinaryQuadraticModel(vartype='BINARY')
    for i in range(N):
        coefficient[i] = -weight[j, i]
    for i in range(M+1):
        coefficient[N+i] = -etaCoefficient[i]
    for i in range(M+1):
        coefficient[N+M+1+j*(M+1)+i] = alphaCoefficient[i, j]

    for i in range(totalNum):
        constraint.add_variable('BINARY', i)
        constraint.set_linear(i, coefficient[i])
    cqm.add_constraint(constraint, sense="==", rhs=b[j])

# volume constraint
constraint = BinaryQuadraticModel(vartype='BINARY')
coefficient = coefficient = np.zeros((totalNum, 1))
for i in range(N):
    coefficient[i] = 1.0 / N
vol = data_dir['targetVol'][0, 0] / N
for i in range(totalNum):
    constraint.add_variable('BINARY', i)
    constraint.set_linear(i, coefficient[i])
cqm.add_constraint(constraint, sense="==", rhs=b[j])

# Need to use your own token here
sampler = LeapHybridCQMSampler(
    token="")
sampleset = sampler.sample_cqm(
    cqm, label='to cqm default')

solution = sampleset.first.sample

x = np.zeros((totalNum, 1))
selected_item_indices = [key for key, val in solution.items() if val == 1.0]
x[selected_item_indices] = 1
x = x[0:N]

sio.savemat('cons_mip_result.mat', {'qaResult': x})
