from VirtualSpinning.Aux import odf 
from matplotlib import pyplot as plt


# Test de Distribucion Random 
def testRand():
   # Create distribution to sample from 
   rd = odf.Rand(lower=-5.) 
   x = [rd() for i in range(1000)]
   xc, dx, h = odf.calcular_histograma(x=x, n=20)
   fig, ax = plt.subplots() 
   ax.bar(xc, h, width=dx, edgecolor='k')
   ax.set_title("testRand()")
   plt.show()


# Test de Distribucion Normal Truncada 
def testNormTr():
   nt = odf.NormTr(loc=1., scale=0.1, lower=1., upper=2.) 
   x = [nt() for i in range(1000)] 
   xc, dx, h = odf.calcular_histograma(x=x, n=20, rango=(1., 2.)) 
   fig, ax = plt.subplots() 
   ax.bar(xc, h, width=dx, edgecolor="k") 
   ax.set_xlim(right=2.)
   ax.set_title("testNormTr()")
   plt.show()


if __name__=="__main__":
    testRand()
    testNormTr()