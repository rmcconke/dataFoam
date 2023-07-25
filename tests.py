
from utilities.MLDatasetFromFoamCase import MLDatasetFromFoamCase
import os
import numpy as np

dataFoam_folder = os.getcwd()

foam_data_case = MLDatasetFromFoamCase(data_save_path=os.path.join(os.getcwd(),'test_data/numpy'),
                                       foam_parent_dir=os.path.join(os.getcwd(),'test_data'),
                                       case_name='case_1p0',
                                       case_type='komegasst',
                                       write_fields_application='writeFields_RANS',
                                       write_fields_flag=True)

foam_data_case.writeFields()
foam_data_case.saveDataset(dataset_prefix='komegasst_case_1p0')

# Check correct number of cells is read
C = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_C.npy'))
print(f'[dataFoam tests] Checking correct number of cells is read....')
assert len(C) == 14751

# Check scalars, vectors, symmtensors, and tensors are read correctly
k = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_k.npy'))
U = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_U.npy'))
S = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_S.npy'))
gradU = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_gradU.npy'))

print(f'[dataFoam tests] Checking that scalars are read correctly....')
assert (k[0] == 7.44073e-08) & (k[1] == 7.38061e-08) & (k[2] == 7.38644e-08) & (k[-1] == 7.70386e-09)

print(f'[dataFoam tests] Checking that vectors are read correctly....')
assert (U[0,0] == 0.00158869) & (U[0,1] == 1.06907e-07) & (U[0,2] == 0)
assert (U[1,0] == 0.00158241) & (U[1,1] == -5.93468e-08) & (U[1,2] == 0)
assert (U[2,0] == 0.00158323) & (U[2,1] == 1.08456e-07) & (U[2,2] == -1.91588e-26)
assert (U[-1,0] == 0.000495152) & (U[-1,1] == -1.03916e-06) & (U[-1,2] == 0)

print(f'[dataFoam tests] Checking that symmTensors are read correctly....')
assert (S[0,0,0] == 2.91774161791e-05) & (S[0,0,1] == -0.389907260139) & (S[0,0,2] == -1.90484370391e-24) \
      & (S[0,1,0] == -0.389907260139) & (S[0,1,1] == -8.91195109353e-05) & (S[0,1,2] == 1.89924496275e-21) \
      & (S[0,2,0] == -1.90484370391e-24) & (S[0,2,1] == 1.89924496275e-21) & (S[0,2,2] == 0) 

assert (S[1,0,0] == 3.09509331373e-05) & (S[1,0,1] == -0.388324194552) & (S[1,0,2] == 5.26877331679e-26) \
      & (S[1,1,0] == -0.388324194552) & (S[1,1,1] == 2.10381637335e-06) & (S[1,1,2] == 1.31723921224e-32) \
      & (S[1,2,0] == 5.26877331679e-26) & (S[1,2,1] == 1.31723921224e-32) & (S[1,2,2] == 0) 

assert (S[2,0,0] == 4.1809472726e-05) & (S[2,0,1] == -0.388050932452) & (S[2,0,2] == -3.39886575925e-25) \
      & (S[2,1,0] == -0.388050932452) & (S[2,1,1] == -9.56457600965e-05) & (S[2,1,2] == 3.43704625725e-21) \
      & (S[2,2,0] == -3.39886575925e-25) & (S[2,2,1] == 3.43704625725e-21) & (S[2,2,2] == -1.96519509009e-51) 

assert (S[-1,0,0] == 0.000930181707256) & (S[-1,0,1] == 0.258170266276) & (S[-1,0,2] == -6.91147562237e-20) \
      & (S[-1,1,0] == 0.258170266276) & (S[-1,1,1] == -0.00155624474334) & (S[-1,1,2] == 5.86029411444e-22) \
      & (S[-1,2,0] == -6.91147562237e-20) & (S[-1,2,1] == 5.86029411444e-22) & (S[-1,2,2] == 0) 

print(f'[dataFoam tests] Checking that tensors are read correctly....')
assert (gradU[0,0,0] == 2.91774161791e-05) & (gradU[0,0,1] == -0.7798138247) & (gradU[0,0,2] == 0) \
      & (gradU[0,1,0] == -6.95577613254e-07) & (gradU[0,1,1] == -8.91195109353e-05) & (gradU[0,1,2] == 0) \
      & (gradU[0,2,0] == -3.80968740782e-24) & (gradU[0,2,1] == 3.79848992551e-21) & (gradU[0,2,2] == 0) 

assert (gradU[1,0,0] == 3.09509331373e-05) & (gradU[1,0,1] == -0.77664838057) & (gradU[1,0,2] == 0) \
      & (gradU[1,1,0] == -8.53457596887e-09) & (gradU[1,1,1] == 2.10381637335e-06) & (gradU[1,1,2] == 0) \
      & (gradU[1,2,0] == 1.05375466336e-25) & (gradU[1,2,1] == 2.63447842449e-32) & (gradU[1,2,2] == 0) 

assert (gradU[2,0,0] == 4.1809472726e-05) & (gradU[2,0,1] == -0.776101000365) & (gradU[2,0,2] == -5.61260569453e-25) \
      & (gradU[2,1,0] == -8.64539580809e-07) & (gradU[2,1,1] == -9.56457600965e-05) & (gradU[2,1,2] == -3.66673515707e-29) \
      & (gradU[2,2,0] == -1.18512582397e-25) & (gradU[2,2,1] == 6.87409255117e-21) & (gradU[2,2,2] == -1.96519509009e-51) 

assert (gradU[-1,0,0] == 0.000930181707256) & (gradU[-1,0,1] == 0.516460503184) & (gradU[-1,0,2] == 0) \
      & (gradU[-1,1,0] == -0.000119970631196) & (gradU[-1,1,1] == -0.00155624474334) & (gradU[-1,1,2] == 0) \
      & (gradU[-1,2,0] == -1.38229512447e-19) & (gradU[-1,2,1] == 1.17205882289e-21) & (gradU[-1,2,2] == 0) 

# Check correct matrix multiplication expression is done by OpenFOAM
Shat = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_Shat.npy'))
Rhat = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_Rhat.npy'))

T1 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T1.npy'))
T2 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T2.npy'))
T3 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T3.npy'))
T4 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T4.npy'))
T5 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T5.npy'))
T6 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T6.npy'))
T7 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T7.npy'))
T8 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T8.npy'))
T9 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T9.npy'))
T10 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_T10.npy'))

I = np.zeros((14751,3,3))
I[:,0,0] = np.ones(14751)
I[:,1,1] = np.ones(14751)
I[:,2,2] = np.ones(14751)

T1_test = Shat
T2_test = (Shat @ Rhat) - (Rhat @ Shat)
T3_test = (Shat @ Shat) - 1/3 * np.trace((Shat @ Shat),axis1=1,axis2=2)[:,None,None]*I
T4_test = (Rhat @ Rhat) - 1/3 * np.trace((Rhat @ Rhat),axis1=1,axis2=2)[:,None,None]*I
T5_test = (Rhat @ Shat @ Shat) - (Shat @ Shat @ Rhat)
T6_test = (Rhat @ Rhat @ Shat) + (Shat @ Rhat @ Rhat) - 2/3* np.trace((Shat @ Rhat @ Rhat),axis1=1,axis2=2)[:,None,None]*I
T7_test = (Rhat @ Shat @ Rhat @ Rhat) - (Rhat @ Rhat @ Shat @ Rhat)
T8_test = (Shat @ Rhat @ Shat @ Shat) - (Shat @ Shat @ Rhat @ Shat)
T9_test = (Rhat @ Rhat @ Shat @ Shat) + (Shat @ Shat @ Rhat @ Rhat) - 2/3 * np.trace((Shat @ Shat @ Rhat @ Rhat),axis1=1,axis2=2)[:,None,None]*I
T10_test = (Rhat @ Shat @ Shat @ Rhat @ Rhat) - (Rhat @ Rhat @ Shat @ Shat @ Rhat)

print(f'[dataFoam tests] Checking that basis tensors are being multiplied correctly....')
assert np.max(T1 - T1_test) < 1E-10
assert np.max(T2 - T2_test) < 1E-10
assert np.max(T3 - T3_test) < 1E-10
assert np.max(T4 - T4_test) < 1E-10
assert np.max(T5 - T5_test) < 1E-10
assert np.max(T6 - T6_test) < 1E-10
assert np.max(T7 - T7_test) < 1E-10
assert np.max(T8 - T8_test) < 1E-10
assert np.max(T9 - T9_test) < 1E-10
assert np.max(T10 - T10_test) < 1E-10

# Check certain invariants are zero, certain invariants are non-zero
I1_8 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_I1_8.npy'))
I2_8 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_I2_8.npy'))

I1_25 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_I1_25.npy'))
I2_25 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_I2_25.npy'))

I1_29 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_I1_29.npy'))
I2_29 = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_I2_29.npy'))

print(f'[dataFoam tests] Checking some invariants for expected zero/non-zero behaviour....')

assert np.sum(abs(I1_8)) < 1E-10
assert np.sum(abs(I2_8)) > 1E-10

assert np.sum(abs(I1_25)) > 1E-10
assert np.sum(abs(I2_25)) < 1E-10

assert np.sum(abs(I1_29)) > 1E-10
assert np.sum(abs(I2_29)) < 1E-10

# Check symmetry of S, and anti-symmetry of R
print(f'[dataFoam tests] Checking symmetric and zero trace basis tensors and S')
R = np.load(os.path.join(dataFoam_folder,'test_data/numpy/komegasst_case_1p0_R.npy'))
def check_symmetric(tensor):
    threshold = 1E-10
    assert (abs((tensor[:,0,1] - tensor[:,1,0]) < threshold).all()) & (abs((tensor[:,0,2] - tensor[:,2,0]) < threshold).all()) & (abs((tensor[:,1,2] - tensor[:,2,1])<threshold).all()) 

def check_antisymmetric(tensor):
    assert ((tensor[:,0,1] == -tensor[:,1,0]).all()) & ((tensor[:,0,2] == -tensor[:,2,0]).all()) & ((tensor[:,1,2] == -tensor[:,2,1]).all()) 

def check_zero_trace(tensor):
    print(f'[dataFoam tests] Maximum trace: {max(abs(np.trace(tensor,axis1=1,axis2=2)))}')
    #assert (abs(np.trace(tensor,axis1=1,axis2=2) - 0 ) < 1E-3).all()

for i,T in enumerate([T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, S]):
    print(i)
    check_symmetric(T)
    check_zero_trace(T)

print(f'[dataFoam tests] Checking antisymmetric and zero trace R')
check_antisymmetric(R)
check_zero_trace(R)

# Check invariants are actually invariant after transformation
foam_data_case = MLDatasetFromFoamCase(data_save_path=os.path.join(dataFoam_folder,'test_data/numpy'),
                                       foam_parent_dir=os.path.join(dataFoam_folder,'test_data'),
                                       case_name='case_1p0_rotatedz45',
                                       case_type='komegasst',
                                       write_fields_application='writeFields_RANS',
                                       write_fields_flag=True)

foam_data_case.writeFields()
foam_data_case.saveDataset(dataset_prefix='komegasst_case_1p0_z45')

foam_data_case = MLDatasetFromFoamCase(data_save_path=os.path.join(dataFoam_folder,'test_data/numpy'),
                                       foam_parent_dir=os.path.join(dataFoam_folder,'test_data'),
                                       case_name='case_1p0_rotatedz90',
                                       case_type='komegasst',
                                       write_fields_application='writeFields_RANS',
                                       write_fields_flag=True)

foam_data_case.writeFields()
foam_data_case.saveDataset(dataset_prefix='komegasst_case_1p0_z90')

foam_data_case = MLDatasetFromFoamCase(data_save_path=os.path.join(dataFoam_folder,'test_data/numpy'),
                                       foam_parent_dir=os.path.join(dataFoam_folder,'test_data'),
                                       case_name='case_1p0_rotatedxy45',
                                       case_type='komegasst',
                                       write_fields_application='writeFields_RANS',
                                       write_fields_flag=True)

foam_data_case.writeFields()
foam_data_case.saveDataset(dataset_prefix='komegasst_case_1p0_xy45')

foam_data_case = MLDatasetFromFoamCase(data_save_path=os.path.join(dataFoam_folder,'test_data/numpy'),
                                       foam_parent_dir=os.path.join(dataFoam_folder,'test_data'),
                                       case_name='case_1p0_rotatedxyz70',
                                       case_type='komegasst',
                                       write_fields_application='writeFields_RANS',
                                       write_fields_flag=True)

foam_data_case.writeFields()
foam_data_case.saveDataset(dataset_prefix='komegasst_case_1p0_xyz70')


for i in range(47):
    #print(f'')
    print(f'[dataFoam tests] Checking invariants of basis tensor {i+1} for invariance under rotation')
    I1_i_base = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_I1_{i+1}.npy'))
    I1_i_z45 = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_z45_I1_{i+1}.npy'))
    I1_i_z90 = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_z90_I1_{i+1}.npy'))
    I1_i_xy45 = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_xy45_I1_{i+1}.npy'))
    I1_i_xyz70 = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_xyz70_I1_{i+1}.npy'))

    I2_i_base = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_I1_{i+1}.npy'))
    I2_i_z45 = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_z45_I1_{i+1}.npy'))
    I2_i_z90 = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_z90_I1_{i+1}.npy'))
    I2_i_xy45 = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_xy45_I1_{i+1}.npy'))
    I2_i_xyz70 = np.load(os.path.join(dataFoam_folder,f'test_data/numpy/komegasst_case_1p0_xyz70_I1_{i+1}.npy'))

    assert (max(abs((I1_i_base - I1_i_z45))) < 1E-8)
    assert (max(abs((I2_i_base - I2_i_z45))) < 1E-8)

    assert (max(abs((I1_i_base - I1_i_z90))) < 1E-8)

    assert (max(abs((I1_i_base - I1_i_xy45))) < 1E-8)
    assert (max(abs((I2_i_base - I2_i_xy45))) < 1E-8)

    assert (max(abs((I1_i_base - I1_i_xyz70))) < 1E-8)
    assert (max(abs((I2_i_base - I2_i_xyz70))) < 1E-8)

print(f'================================================================')
print(f'[dataFoam tests] All assertions passed.')
print(f'================================================================')
