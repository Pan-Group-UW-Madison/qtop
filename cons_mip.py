import gurobipy as gp
from gurobipy import GRB
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt


def populate_and_solve(m):
    data_dir = sio.loadmat('cons_mip.mat')
    weight = np.array(data_dir['weight'])
    b = np.array(data_dir['b'])
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
    print(Nc)

    objExpr = eta

    A = 2 * Upper
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
        subConsExpr += rhoTuple[i]
    subConsExpr -= data_dir['vol'][0, 0]
    consExpr += A * (subConsExpr)**2

    m.setObjective(objExpr + consExpr, GRB.MINIMIZE)
    m.Params.MIPGap = 1e-3
    m.optimize()

    eta = Upper*(1 - 1 / (2**N))*etaTuple[0].X
    for i in range(1, M+1):
        eta += Upper/(2**i)*etaTuple[i].X

    x = np.zeros((N, 1))
    for i in range(0, N):
        x[i] = rhoTuple[i].X
    sio.savemat('cons_mip_result.mat', {'x': x})

    # fig = plt.figure()
    # fig.set_figheight(2)
    # fig.set_figwidth(6)
    # plt.imshow(1 - x.reshape(60, 20).T, cmap='gray', vmin=0, vmax=1)
    # plt.axis("off")
    # fig.tight_layout()
    # plt.savefig("cons_mip.eps")


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
        populate_and_solve(model)
