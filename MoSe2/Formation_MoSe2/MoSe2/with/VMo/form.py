import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-5997075.714803956
E_pristine=-6108318.476328624
mu_Mo=-222473.348568306/2
Ef=E_defect-E_pristine+(1*mu_Mo)
print(Ef)
