import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import pickle
import pandas as pd

import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'

def classify(model_path, rdf):
    print("start!")
    model = DSFNet()
    print("passed  model = DSFNet()")
    model.load_state_dict(torch.load(model_path))
    print("passed  model.load...")
    model.eval()
    loader = torch.utils.data.DataLoader(format_data(rdf), batch_size=512, shuffle=False)
    predicted = []
    wells = []
    for data in loader:
      predicted.append(model(data[0]).cpu().detach())
      wells.append(data[2])

    predicted = np.hstack(predicted)
    wells = np.hstack(wells)
    df = pd.DataFrame(predicted, index=wells, columns=['Prob 0', 'Prob 1', 'Prob 2'])
    df['well'] = df.index
    print(df)
    return df

def format_data(rdf):
    data_array = [np.array(rdf[k]) for k in sorted(rdf.keys()) if k != 'well']
    data_array = np.vstack(data_array)
    wells = np.array(rdf['well'])
    data_array = data_array.T
    data_array = data_array.astype('float32')
    data_array = normalize_by_row(data_array)
    data_array = np.array([row.reshape(12, 70) for row in data_array])
    print(data_array.shape)
    y = np.zeros((data_array.shape[0]))
    return ArrayDataset(data_array, y, wells)
	

def normalize_by_row(matrix):
    max_of_rows = np.max(matrix, axis=1)
    return matrix/max_of_rows[:,None]
    

class ArrayDataset(torch.utils.data.Dataset):
    """
    Operates on Numpy Arrays for Transforms
    """
    def __init__(self,
                 X_array,
                 y_array,
                 names,
                 transform=None):
        assert X_array.shape[0] == y_array.shape[0]
        self.X = X_array
        self.y  = y_array
        self.names = names
        self.transform = transform
        
    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        X_item = self.X[idx]
        y_item = self.y[idx]
        name_item = self.names[idx]
        if self.transform:
            sample = self.transform(X_item)
        return X_item, y_item, name_item


class DSFNet(nn.Module):
    """
    Doc string what?
    """

    def __init__(self,
                 input_dim=(12, 70)):
        super().__init__()

        self.conv1 = nn.Conv1d(12, 32, 3, stride=1, padding=0)  # 12 x 70 --> 32 x 68
        self.conv2 = nn.Conv1d(32, 64, 3, stride=1, padding=0)  # 32 x 34 --> 64 x 32 (mp to 16)

        self.fc1 = nn.Linear(64 * 16, 50)
        # self.fc2 = nn.Linear(50, 20)
        self.out = nn.Linear(50, 3)

    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = F.max_pool1d(x, 2)
        x = F.relu(self.conv2(x))
        x = F.max_pool1d(x, 2)
        # x = F.relu(self.conv3(x))
        # x = F.max_pool1d(x, 2)
        x = x.view(-1, 64 * 16)
        x = F.relu(self.fc1(x))
        x = self.out(x)
        return F.log_softmax(x, dim=1)


def write_data(rdf):
  with open('toy_data.pickle', 'wb') as handle:
    pickle.dump(rdf, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
def load_data():
  with open('toy_data.pickle', 'rb') as handle:
   return pickle.load(handle)



