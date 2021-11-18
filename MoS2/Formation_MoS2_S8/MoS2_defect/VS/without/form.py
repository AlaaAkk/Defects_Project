import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-3314031.328428974
E_pristine=-3324912.918890152

mu=-87031.346957329/8


Ef=E_defect-E_pristine+(1*mu)
print(Ef)
