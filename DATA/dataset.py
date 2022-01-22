from __future__ import print_function
import torch.utils.data as data
import scipy.io
import numpy as np

def data_load(data_name):
    data_dir = "./data/BIC/"
    # data_name = "SARCOMA_Gene_Expression.txt"
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
class EYB(data.Dataset):
    def __init__(self, transform=None):
        self.transform = transform

        self.train_num = int(624)

        data_dir = "./data/BIC/"
        data_name1 = "BREAST_Gene_Expression.txt"
        data_name2 = "BREAST_Methy_Expression.txt"
        data_name3 = "BREAST_Mirna_Expression.txt"
        gene_expression1, Img1, Label1, n_input1 = load_data(data_name1)
        gene_expression2, Img2, Label2, n_input2 = load_data(data_name2)
        gene_expression3, Img3, Label3, n_input3 = load_data(data_name3)



        Img3 = Img3[:, :, 0:885, :]
        Img1 = np.squeeze(Img1, axis=None)
        Img2 = np.squeeze(Img2, axis=None)
        Img3 = np.squeeze(Img3, axis=None)
        self.Img1 = Img1.astype(np.float32)
        self.Img2 = Img2.astype(np.float32)
        self.Img3 = Img3.astype(np.float32)
        print(self.Img1.shape)
        print(self.Img2.shape)
        print(self.Img3.shape)

    def __getitem__(self, index):

        img_traina, img_trainb, img_trainc = self.Img1[index, :], self.Img2[index, :], self.Img3[index, :]
        return img_traina, img_trainb,img_trainc

    def __len__(self):
        return self.train_num



