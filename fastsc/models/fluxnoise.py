import numpy as np


def get_flux_noise(omega, device):
    # Return the deviation of frequency due to gaussian noise in flux
    # omega is the target frequency, device.flux_sigma is the std dev.
    Ejs = device.ejs
    Ejl = device.ejl
    Ec = device.ec
    gamma1 = Ejl/Ejs
    d1 = (gamma1 - 1)/(gamma1 + 1)
    high = float(device.omega_max) + 1e-5
    low = float(device.omega_min) - 1e-5
    print(omega, high, low)
    if (omega < low or omega > high):
        print("Flux noise warning: " + str(omega) + " not in range [" + str(low) + "," + str(high) + "].")
        return omega
    flux_est = np.arccos((omega-(high+low)/2.0)/((high-low)/2.0))/(2*np.pi)
    flux_new = flux_est + np.random.normal(0,device.flux_sigma)
    def _Em(phi, m):
        # Asymmetric transmon
        return -(Ejs + Ejl)/2 + np.sqrt(4*(Ejs + Ejl)*Ec*np.sqrt(np.cos(np.pi*phi)**2 + d1**2*np.sin(np.pi*phi)**2))*(m + 1/2) - Ec*(6*m**2 + 6*m + 3)/12
    omega_old = _Em(flux_est, 1) - _Em(flux_est, 0)
    omega_new = _Em(flux_new, 1) - _Em(flux_new, 0)
    if omega_new > high:
        omega_new = high
    elif omega_new < low:
        omega_new = low
    return omega_new - omega_old
