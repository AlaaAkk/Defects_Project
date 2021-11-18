import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-14192154.320500717
E_pristine=-14644858.463124419
mu_W=-905397.333160509/2
Ef=E_defect-E_pristine+(1*mu_W)
print(Ef)
