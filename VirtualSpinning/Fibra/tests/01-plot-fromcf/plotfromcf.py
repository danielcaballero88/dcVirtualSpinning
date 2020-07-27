from VirtualSpinning.Fibra.Fibra import Fibra
from matplotlib import pyplot as plt
import os.path


CWD = os.path.dirname(__file__)


def main():
    cfile = os.path.join(CWD, "ConfigurationFile.txt")
    #
    fib = Fibra.from_cf(cfile, lamr=1.0)
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