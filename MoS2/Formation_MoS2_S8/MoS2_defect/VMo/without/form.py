import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-3213669.749963366
E_pristine=-3324912.918890152

mu=-222471.853000009/2


Ef=E_defect-E_pristine+(1*mu)
print(Ef)
