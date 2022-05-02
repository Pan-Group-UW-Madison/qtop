import gurobipy as gp
from gurobipy import GRB
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt


def populate_and_solve(m):
    data_dir = sio.loadmat('mip.mat')
    f = data_dir['f']

    N = len(f[0])
    rhoTuple = m.addVars(N, vtype=GRB.BINARY)

    objExpr = 0
    for i in range(N):
        objExpr += f[0, i] * rhoTuple[i]

    A = 2 * data_dir['obj'][0, 0]
    consExpr = 0
    for i in range(0, N):
        consExpr += rhoTuple[i]
    consExpr -= data_dir['vol'][0, 0]
    consExpr = A * (consExpr)**2

    m.setObjective(objExpr + consExpr, GRB.MINIMIZE)
    m.Params.MIPGap = 1e-2
    m.optimize()

    x = np.zeros((N, 1))
    for i in range(0, N):
        x[i] = rhoTuple[i].X
    sio.savemat('mip_result.mat', {'x': x})

    # fig = plt.figure()
    # fig.set_figheight(2)
    # fig.set_figwidth(6)
    # plt.imshow(1 - x.reshape(60, 20).T, cmap='gray', vmin=0, vmax=1)
    # plt.axis("off")
    # fig.tight_layout()
    # plt.savefig("lin_mip.eps")


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
