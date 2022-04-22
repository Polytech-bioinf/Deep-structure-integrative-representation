import torch.nn as nn
import torch


class Networks(nn.Module):
    def __init__(self):
        super(Networks, self).__init__()


        self.encoder1 = nn.Sequential(
            nn.Conv2d(1, 15, kernel_size=(5, 1), stride=5, padding=0, bias=True),
            nn.ReLU(),

        )
        self.decoder1 = nn.Sequential(
            nn.ConvTranspose2d(15, 1, kernel_size=(5, 1), stride=5, padding=0, bias=False),
            nn.ReLU(),

        )

        self.encoder2 = nn.Sequential(
            nn.Conv2d(1, 15, kernel_size=(5, 1), stride=5, padding=0, bias=True),
            nn.ReLU(),

        )
        self.decoder2 = nn.Sequential(

            nn.ConvTranspose2d(15, 1, kernel_size=(5, 1), stride=5, padding=0, bias=False),
            nn.ReLU(),
        )
        self.encoder3 = nn.Sequential(
            nn.Conv2d(1, 15, kernel_size=(5, 1), stride=5, padding=0, bias=False),
            nn.ReLU(),

        )
        self.decoder3 = nn.Sequential(

            nn.ConvTranspose2d(15, 1, kernel_size=(5, 1), stride=5, padding=0, bias=False),
            nn.ReLU(),
        )

        self.weight = nn.Parameter(1.0e-4 * torch.ones(624, 624))




    def forward(self, input1, input2,input3):
        output1 = self.encoder1(input1)
        output1 = self.decoder1(output1)

        output2 = self.encoder2(input2)
        output2 = self.decoder2(output2)

        output3 = self.encoder3(input3)
        output3 = self.decoder3(output3)
        return output1, output2, output3

    def forward2(self, input1, input2, input3):
        #coef = abs(self.weight - torch.diag(torch.diag(self.weight)))
        coef = (self.weight - torch.diag(torch.diag(self.weight)))
        z1 = self.encoder1(input1)
        z1 = z1.view(624, -1)
        zcoef1 = torch.matmul(coef, z1)
        output1 = zcoef1.view(624, 15, -1, 1)
        output1 = self.decoder1(output1)

        z2 = self.encoder2(input2)
        z2 = z2.view(624, -1)
        zcoef2 = torch.matmul(coef, z2)
        output2 = zcoef2.view(624, 15, -1, 1)
        output2 = self.decoder2(output2)

        z3 = self.encoder3(input3)
        z3 = z3.view(624, -1)
        zcoef3 = torch.matmul(coef, z3)
        output3 = zcoef3.view(624, 15, -1, 1)
        output3 = self.decoder3(output3)
        return z1, z2, z3, output1, output2, output3, zcoef1, zcoef2, zcoef3, coef






