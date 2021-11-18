import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-14578289.161039500
E_pristine=-14644838.711835997
mu_Se=-532375.334009259/8
Ef=E_defect-E_pristine+(1*mu_Se)
print(Ef)
