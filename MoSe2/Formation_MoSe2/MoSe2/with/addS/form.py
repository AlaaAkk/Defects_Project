import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-6174864.306138340
E_pristine=-6108318.476328624

mu=-532375.657399940/8


Ef=E_defect-E_pristine-(1*mu)
print(Ef)
