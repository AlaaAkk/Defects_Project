import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-11850577.239372833
E_pristine=-11861458.807127204

mu=-87031.346957329/8


Ef=E_defect-E_pristine+(1*mu)
print(Ef)
