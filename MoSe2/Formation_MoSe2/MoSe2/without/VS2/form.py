import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-5975198.692503990
E_pristine=-6108297.545624177
mu_Se=-532375.334009259/8
Ef=E_defect-E_pristine+(2*mu_Se)
print(Ef)
