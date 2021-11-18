import matplotlib.pyplot as plt
from numpy import array, reshape, zeros, append, arange
E_defect=-5997055.608034546
E_pristine=-6108297.545624177

mu=-222471.853000009/2


Ef=E_defect-E_pristine+(1*mu)
print(Ef)
