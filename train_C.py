import torch
from data.dataloader import data_loader_train
from models.network import Networks
import models.metrics as metrics
import numpy as np
import scipy.io as sio
from scipy.sparse.linalg import svds
from sklearn import cluster
from sklearn.preprocessing import normalize
import scipy.io
import os
import tensorflow as tf

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")



def data_load(data_name):
    
    data_dir = "./data/BIC/"
    # data_name = "BIC_Gene_Expression.txt"
    data_path = data_dir + data_name
    f = open(data_path)  # 读取文件
    lines = f.readlines()
    # print(len(lines))

    for i in range(len(lines)):
        lines[i] = lines[i].strip().split("\t")[1:]
        # print(lines[i])

    gene_expression = np.zeros([len(lines) - 1, len(lines[1])])

    for i in range(len(lines)):
        if i == 0:
            continue
        for j in range(len(lines[1])):
            gene_expression[i - 1][j] = float(lines[i][j])

    gene_expression = np.matrix(gene_expression)
    # print(gene_expression.shape)
    gene_expression = np.transpose(gene_expression)
    gene_expression = np.array(gene_expression)

    return gene_expression


def load_data(data_name):
    gene_expression = data_load(data_name)
    Label = gene_expression
    # Img = np.reshape(Img, (Img.shape[0], 32, 32, 1))
    Img = np.reshape(gene_expression, [gene_expression.shape[0], 1, gene_expression.shape[1], 1])
    n_input = [1, gene_expression.shape[1]]

    return gene_expression, Img, Label, n_input


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
tf.logging.set_verbosity(tf.logging.ERROR)

data_dir = "./data/BIC/"
data_name1 = "BREAST_Gene_Expression.txt"
data_name2 = "BREAST_Methy_Expression.txt"
data_name3 = "BREAST_Mirna_Expression.txt"
gene_expression1, Img1, Label1, n_input1 = load_data(data_name1)
gene_expression2, Img2, Label2, n_input2 = load_data(data_name2)
gene_expression3, Img3, Label3, n_input3 = load_data(data_name3)
batch_size = gene_expression1.shape[0]
Img3 = Img3[:, :, 0:885, :]
Img1 = np.squeeze(Img1, axis=None)
Img2 = np.squeeze(Img2, axis=None)
Img3 = np.squeeze(Img3, axis=None)
from bunch import *
s = Bunch()
s.data = [[Img1], [Img2], [Img3]]


regg2 = [1, 5, 10, 20, 50, 100, 150, 200]
regg3 = [1, 5, 10, 20, 50, 100, 150, 200]
model = Networks().to(device)
criterion = torch.nn.MSELoss(reduction='sum')
optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=0.001, weight_decay=0.0)
n_epochs = 50
for epoch in range(n_epochs):
    for data in data_loader_train:
        train_imga, train_imgb, train_imgc = data

        input1 = train_imga.view(624, 1, 2000, 1).to(device)
        input2 = train_imgb.view(624, 1, 2000, 1).to(device)
        input3 = train_imgc.view(624, 1, 885, 1).to(device)

        output1, output2, output3 = model(input1, input2, input3)
        loss = criterion(output1, input1) + criterion(output2, input2) + criterion(output3, input3)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    if epoch % 10 == 0:
        print("Epoch {}/{}".format(epoch, n_epochs))
        print("Loss is:{:.4f}".format(loss.item()))
torch.save(model.state_dict(), './models/BIC.pth')

import snf
affinity_networks = snf.make_affinity(s.data, metric='euclidean', K=20, mu=0.5)
fused_network = snf.snf(affinity_networks, K=20)
print(fused_network)
n_epochs2 = 200
optimizer2 = torch.optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=0.001, weight_decay=0.0)
criterion3 = torch.nn.MSELoss(reduction='sum')
criterion2 = torch.nn.L1Loss(reduction='sum')

y = 0  # 记录reg组合
for r in range(0, len(regg2)):
    reg2 = regg2[r]
    t = 0
    for h in range(0, len(regg3)):
        model.load_state_dict(torch.load('./models/BIC.pth', map_location=torch.device('cpu')))
        optimizer2 = torch.optim.Adam(filter(lambda p: p.requires_grad, model.parameters()), lr=0.001, weight_decay=0.0)
        print("step2")
        print("---------------------------------------")
        reg3 = regg3[h]
        l = 0
        for epoch in range(n_epochs2):
            for data in data_loader_train:

                train_imga, train_imgb, train_imgc = data

                input1 = train_imga.view(624, 1, 2000, 1).to(device)
                input2 = train_imgb.view(624, 1, 2000, 1).to(device)
                input3 = train_imgc.view(624, 1, 885, 1).to(device)

                z1, z2, z3, output1, output2, output3, zcoef1, zcoef2, zcoef3, coef = model.forward2(input1, input2,
                                                                                                     input3)
                loss_re = criterion3(coef, torch.zeros(624, 624, requires_grad=True).to(device))
                loss_e = criterion3(zcoef2, z2) + criterion3(zcoef1, z1) + criterion3(zcoef3, z3)
                loss_r = criterion3(output1, input1) + criterion3(output2, input2)+ criterion3(output3, input3)
                loss_ff = criterion2(coef.mul(torch.tensor(fused_network)), torch.zeros((624, 624)))
                loss = loss_r + 1 * loss_e + reg2 * loss_re - reg3 * loss_ff

                optimizer2.zero_grad()
                loss.backward()
                optimizer2.step()
            if epoch % 10 == 0:

                print("Epoch {}/{}".format(epoch, n_epochs2))
                print("Loss is:{:.4f}".format(loss.item()))
                commonC = coef.cpu().detach().numpy()


                np.savetxt("/Users/yangyan/Desktop/test/BIC/_param_comb/param2_" + str(regg2[y]) + "/param3_" + str(regg3[t]) + "/epoch_" + str(l) + "/matrix" + ".txt", np.c_[commonC],
                           fmt='%lf', delimiter='\t')
                torch.save(model.state_dict(), './models/BIC/param2_'+ str(regg2[y]) + '/param3_' + str(regg3[t]) + '/epoch_' + str(l*10) +'.pth')
                l = l + 1
        t = t + 1
    y = y + 1


