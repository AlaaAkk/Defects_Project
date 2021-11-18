import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-3213688.902122602
E_pristine=-3324932.942491689

mu=-222473.348568306/2


Ef=E_defect-E_pristine+(1*mu)
print(Ef)
