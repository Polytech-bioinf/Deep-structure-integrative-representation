import torch
from data.dataset import EYB
data_loader_train = torch.utils.data.DataLoader(EYB(), batch_size=624, shuffle=False)



