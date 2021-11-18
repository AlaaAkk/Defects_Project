import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-11872355.158936081
E_pristine=-11861477.311099075

mu=-87031.629750345/8


Ef=E_defect-E_pristine-(1*mu)
print(Ef)
