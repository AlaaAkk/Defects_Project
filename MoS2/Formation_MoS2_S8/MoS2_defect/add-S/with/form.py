import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-3335811.376811718
E_pristine=-3324932.942491689

mu=-87031.629750345/8


Ef=E_defect-E_pristine-(1*mu)
print(Ef)
