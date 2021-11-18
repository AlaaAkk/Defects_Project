import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-11839695.585247073
E_pristine=-11861458.807127204

mu=-87031.346957329/8


Ef=E_defect-E_pristine+(2*mu)
print(Ef)
