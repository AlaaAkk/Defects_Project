import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-5975218.839354995
E_pristine=-6108318.476328624

mu=-532375.657399940/8


Ef=E_defect-E_pristine+(2*mu)
print(Ef)
