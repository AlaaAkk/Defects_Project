import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-3303169.071443753
E_pristine=-3324932.942491689

mu=-87031.629750345/8


Ef=E_defect-E_pristine+(2*mu)
print(Ef)
