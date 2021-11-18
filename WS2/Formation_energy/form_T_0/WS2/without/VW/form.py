import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-11408754.253221687
E_pristine=-11861458.807127204

mu=-905396.073585909/2


Ef=E_defect-E_pristine+(1*mu)
print(Ef)
