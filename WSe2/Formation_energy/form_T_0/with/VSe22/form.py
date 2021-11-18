import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-14511758.944663303
E_pristine=-14644858.463124419
mu_Se=-532375.657399940/8
Ef=E_defect-E_pristine+(2*mu_Se)
print(Ef)
