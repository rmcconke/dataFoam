import numpy as np
import os

def load_dataset_fields(dataset_folder,genmethod,cases,field):
    data = np.concatenate([np.load(os.path.join(
        dataset_folder,genmethod,genmethod+'_'+case+'_'+field + '.npy')) for case in cases])
    return data


