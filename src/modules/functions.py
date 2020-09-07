from modules.constants import *
#
def f2F(t,f):
    F = 0.0*f
    dt = t[1] - t[0]
    Nt = len(t)
    for it in range(Nt-1):
        F[it + 1] = F[it] + (f[it] + f[it + 1])*dt/2.0
    F[Nt - 1] = 2.0*F[Nt - 2] - F[Nt - 3]
    return F
def F2f(t,F):
    f = 0.0*F
    dt = t[1] - t[0]
    Nt = len(t)
    for it in range(1,Nt-1):
        f[it] = (F[it + 1] - F[it - 1])/(2.0*dt)
    f[0] = 2.0*f[1] - f[2]
    f[Nt - 1] = 2.0*f[Nt - 2] - f[Nt - 3]
    return f
#
def Make_fields(PC):
    T = PC.Ncycle*(tpi/PC.omegac)
    t = np.linspace(0.0, T, PC.Nt)
    A = (PC.E0/PC.omegac)*np.sin(PC.omegac*t)
    E = -PC.E0*np.cos(PC.omegac*t)
    return t, A, E
#
def k2ec(PC,k):
    ec = 0.5*k**2/PC.mc + PC.Eg
    return ec
def k2ev(PC,k):
    ev = 0.5*k**2/PC.mv
    return ev
def k2vc(PC,k):
    vc = k/PC.mc
    return vc
def k2vv(PC,k):
    vv = k/PC.mv
    return vv
#
def get_kpAt(PC,A):
    kpAt = np.zeros([PC.Nk,PC.Nt], dtype='float64')
    for ik in range(PC.Nk):
        kpAt[ik, :] = PC.k[ik] + A[:]
    return kpAt
def ekt2thetakt(PC, t, ekt):
    thetakt = 0.0*ekt
    for ik in range(PC.Nk):
        thetakt[ik, :] = f2F(t,ekt[ik, :])
    return thetakt
#
def get_dcckt(PC, thetavkt, thetackt, evkt, eckt, E):
    dcckt = (0.0 + 0.0*zI)*thetackt
    for ik in range(PC.Nk):
        dcckt[ik, :] = np.exp(zI*(thetackt[ik, :] - thetavkt[ik, :]))*(-E[:])/(eckt[ik, :] - evkt[ik, :])
    return dcckt
#
def dcckt2cckt(PC, t, dcckt):
    cckt = 0.0*dcckt
    for  ik in range(PC.Nk):
        cckt[ik, :] = f2F(t,dcckt[ik, :])
    return cckt
def get_J(PC, cckt,thetavkt,thetackt):
    jk = 0.0*cckt
    for ik in range(PC.Nk):
        jk[ik, :] = cckt[ik, :]* np.exp(zI*(thetavkt[ik, :] - thetackt[ik, :]))
    J = np.sum(jk, axis=0)/PC.Nk
    J = 2.0*np.real(J)
    return J
#
#
def k02k0pAt(PC, k0, A, Ndim=1):
    if (Ndim == 1):
        k0pAt = k0 + A[:]
    elif (Ndim == 2):
        k0pAt = np.zeros([2,PC.Nt],dtype='float64')
        k0pAt[0,:] = k0[0] + 1.0*A[:]
        k0pAt[1,:] = k0[1] + 0.0*A[:]
    return k0pAt

def et2thetat(PC, t, et):
    thetat = f2F(t,et[:])
    return thetat
#
def get_dcct(PC, thetavt, thetact, evt, ect, E):
    dcct = np.exp(zI*(thetact[:] - thetavt[:]))*(-E[:])/(ect[:] - evt[:])
    return dcct
#
def dcct2cct(PC, t, dcct):
    cct = f2F(t,dcct[:])
    return cct
#