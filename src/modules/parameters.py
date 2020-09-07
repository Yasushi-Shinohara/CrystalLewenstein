from modules.constants import *

class parameter_class:
    def __init__(self):
        self.Nk = 20
        self.Nk2D = None
        self.Nk3D = None
        self.a = 10.0
        self.b = None
        self.k = None
        self.k2D = None
        self.k3D = None
        self.mv = -0.5
        self.mc = 0.1
        self.Eg = 0.2
        self.Nt = 1000
        self.Ncycle = 4
        self.E0 = 0.01
        self.omegac = 1.55/Hartree

    def Make_kspace(self, Ndim=1):
        self.b = tpi/self.a
        self.k = np.linspace(-0.5*self.b, 0.5*self.b, self.Nk)
        if (Ndim == 2):
            k = np.linspace(-0.5*self.b, 0.5*self.b, self.Nk)
            self.Nk2D = self.Nk**2
            self.k2D = np.zeros([self.Nk2D, 2], dtype='float64')
            nk = 0
            for ik in range(self.Nk):
                for jk in range(self.Nk):
                    nk = ik + jk*self.Nk
                    self.k2D[nk, 0] = self.k[ik]
                    self.k2D[nk, 1] = self.k[jk]
                    #print(self.k2D[nk,:]/self.b)
                    #print(ik, jk, nk, self.Nk**2)
        #elif (Ndim == 3):
