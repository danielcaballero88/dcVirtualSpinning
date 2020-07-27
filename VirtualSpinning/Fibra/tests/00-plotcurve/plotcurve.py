from VirtualSpinning.Fibra.Fibra import Fibra
from matplotlib import pyplot as plt


def main():
    Et = 3.0e3 
    EbEt = 1.0e-3 
    doteps = 1.0e-3 
    s0 = 17.0e0
    nh = 0.5
    tenbrk = 100.
    lamr = 1. 
    #
    fib = Fibra(Et, EbEt, doteps, s0, nh, lamr, tenbrk)
    #
    out = fib.traccionar(1.0e-2, 3.5e-3, 1.2)
    # 
    fig, ax = plt.subplots(figsize=(12,9))
    # 
    ax.plot(out['time'], out['ten'])
    #
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()