from pycbc import types, fft, waveform
from astropy.cosmology import Planck18  # Planck 2018
c = 3e8
G = 6.67e-11
pc2m = 3.086e16
yr2s = 3.154e7
solarmass = 2e30

def twomass_localrate(m1,m2,p_m1):
    # fix the primary mass's f_pbh to be 3e-3
    f_pbh_m1 = 3e-3
    # the fraction of p1 is a free variable
    p1 = p_m1
    p2 = 1-p1
    # the total fraction of PBH including m1 and m2
    f = f_pbh_m1/p_m1
    sigma_eq = 0.005
    return 3e6*f**2*(0.7*f**2+sigma_eq**2)**(-21./74) \
            *min(p1/m1,p2/m2)*(p1/m1+p2/m2)*(m1*m2)**(3./37)*(m1+m2)**(36./37)


def dE_df(mass1,mass2,z,sens_freq=40):
    '''
    Return the dE/df given a mass and distance of a BBH merger
    '''
    distance_mpc = Planck18.luminosity_distance(z).value
    f_source = (1+z) * sens_freq
    if f_source > 2048:
        return 0
    else:
        delta_f = 1/4
        sptilde, _ = waveform.get_fd_waveform(approximant="IMRPhenomD",
                                              mass1=mass1, 
                                              mass2=mass2, 
                                              delta_f=1/4,
                                              f_higher = 2048,
                                              f_lower=sens_freq,
                                              distance=distance_mpc)
        kind = int(f_source / delta_f)
        try:
            sptildedata = sptilde.data[kind]
        except IndexError:
            return 0
        return 16/5*np.pi**2*c**3/G*f_source**2*(distance_mpc*1e6*pc2m)**2* \
                (np.abs(sptildedata))**2

t0 = Planck18.age(0).value

def twomass_rate(m1,m2,p_m1,z):
    f_pbh_m1 = 3e-3
    p1 = p_m1
    p2 = 1-p1
    f = f_pbh_m1/p_m1
    sigma_eq = 0.005
    t = Planck18.age(z).value
    return 3e6 * (t/t0)**(-34/37) * f**2*(0.7*f**2+sigma_eq**2)**(-21./74)   \
            *min(p1/m1,p2/m2)*(p1/m1+p2/m2)    \
            *(m1*m2)**(3./37)*(m1+m2)**(36./37)

def E_of_z(z):
    return np.sqrt(Planck18.Om0*(1+z)**3+(1-Planck18.Om0))


cont_rate = (1e9*pc2m)**3*yr2s
def int_of_z(z,mass1,mass2,p_m1,sens_freq=40):
    return twomass_rate(mass1,mass2,p_m1,z)/cont_rate/(1+z)/E_of_z(z)*dE_df(mass1,mass2,z,sens_freq)

rho_c = 3*lal.H0_SI**2*c**2/8/np.pi/G