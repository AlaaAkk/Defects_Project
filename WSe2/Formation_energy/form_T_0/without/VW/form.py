import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-14192135.482328624
E_pristine=-14644838.711835997
mu_Mo=-905396.073585909/2
Ef=E_defect-E_pristine+(1*mu_Mo)
print(Ef)
