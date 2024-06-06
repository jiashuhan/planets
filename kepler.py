import numpy as np

Msun = 1.9884e30 # [kg]
G    = 6.6743e-11 # [m^3/kg/s^2]

def kep2cart(e, a, i, Ω, m, ϖ=None, ω=None, L=None, M=None, m0=Msun):
    """
    Convert Keplerian orbital elements of a given epoch to state vectors in 
    the heliocentric cartesian coordinates in the ecliptic plane. The x-axis
    points toward the vernal equinox (in the solar system case).

    Parameters
    ----------
    e: float
        Eccentricity [rad]
    a: float
        Semi-major axis [m]
    i: float
        Inclination [deg]
    Ω: float
        Longitude of the ascending node [deg]
    m: float
        Mass of orbiting object [kg]
    ϖ: float
        Longitude of periapsis [deg]; used if ω is None
    ω: float
        Argument of periapsis [deg]
    L: float
        Mean longitude [deg] at given epoch; used if M is None
    M: float
        Mean anomaly [deg] at given epoch
    m0: float, optional
        Mass of central body [kg], set to 1 solar mass by default

    Returns
    -------
    x, y, z, vx, vy, vz: tuple
        Position [m] and velocity [m/s] measured in the 
        initial rest frame of the central body.

    Notes
    -----
    (1) e and a define the shape and size of the ellipse. i and 
    Ω define the orientation of the orbital plane. ω or ϖ defines
    the orientation of the ellipse in the orbital plane and 
    M or L defines the position of the orbiting body along the
    ellipse at a specific time. Alternatively, one can use q
    (periapsis) in place of a, or T (time of perihelion passage) 
    in place of M or L.

    (2) For solar system planets, the elements used are e, a, i, Ω, ϖ, L.
    """
    if ω is None: # provided ϖ
        ω = (ϖ - Ω) * np.pi / 180 # argument of periapsis [rad]
    else: # provided ω
        ω = ω * np.pi / 180
    
    i = i * np.pi / 180
    Ω = Ω * np.pi / 180

    if M is None:
        M = L - ϖ # mean anomaly [deg]

    E  = ecc_anomaly(e, M) * np.pi / 180 # eccentric anomaly [rad]
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2)) # true anomaly [rad]

    # heliocentric coordinates [m] in orbital plane, z0 = 0
    # r0 = a * (1 - e * np.cos(E))
    x0 = a * (np.cos(E) - e) # = r0 * np.cos(nu)
    y0 = a * (1 - e**2)**0.5 * np.sin(E) # = r0 * np.sin(nu)

    # ecliptic coordinates [m]
    x =   ( np.cos(ω) * np.cos(Ω) - np.sin(ω) * np.sin(Ω) * np.cos(i)) * x0 \
        + (-np.sin(ω) * np.cos(Ω) - np.cos(ω) * np.sin(Ω) * np.cos(i)) * y0
    y =   ( np.cos(ω) * np.sin(Ω) + np.sin(ω) * np.cos(Ω) * np.cos(i)) * x0 \
        + (-np.sin(ω) * np.sin(Ω) + np.cos(ω) * np.cos(Ω) * np.cos(i)) * y0
    z =   np.sin(ω) * np.sin(i) * x0 + np.cos(ω) * np.sin(i) * y0

    # find velocity in ecliptic coordinates [m/s]
    mu = G * (m + m0)
    v  = np.sqrt(mu / a * (1 + e * np.cos(E)) / (1 - e * np.cos(E))) # orbital speed [m/s]
    
    theta = np.arcsin(np.sqrt((1 - e**2) / (1 - e**2 * np.cos(E)**2))) # angle between r and v [rad]
    vx0, vy0 = v * np.cos(nu + theta), v * np.sin(nu + theta) # velocity in orbital plane [m/s]

    vx =   ( np.cos(ω) * np.cos(Ω) - np.sin(ω) * np.sin(Ω) * np.cos(i)) * vx0 \
         + (-np.sin(ω) * np.cos(Ω) - np.cos(ω) * np.sin(Ω) * np.cos(i)) * vy0
    vy =   ( np.cos(ω) * np.sin(Ω) + np.sin(ω) * np.cos(Ω) * np.cos(i)) * vx0 \
         + (-np.sin(ω) * np.sin(Ω) + np.cos(ω) * np.cos(Ω) * np.cos(i)) * vy0
    vz =   np.sin(ω) * np.sin(i) * vx0 + np.cos(ω) * np.sin(i) * vy0

    return (x, y, z), (vx, vy, vz)

# Solve for eccentric anomaly E [deg] given eccentricity e [rad] and mean anomaly M [deg] using Newton's method
def ecc_anomaly(e, M, tolerance=0.002):
    e_deg = e * 180 / np.pi # e [deg]
    E  = M + e_deg * np.sin(M * np.pi / 180) # initial guess for E [deg]
    dE = 1

    while np.abs(dE/E) > tolerance:
        dM = M - (E - e_deg * np.sin(E * np.pi / 180))
        dE = dM / (1 - e * np.cos(E * np.pi / 180))
        E += dE

    return E # [deg]

# Mean anomaly [deg] given eccentricity [rad] and eccentric anomaly [deg]
def mean_anomaly(e, E):
    return E - e * 180 / np.pi * np.sin(E * np.pi / 180)

# Converts state vectors (r [m] and v [m/s]) into orbital elements
def cart2kep(x, y, z, vx, vy, vz, m, m0=Msun):
    mu  = G * (m + m0) # [m^3/s^2]
    r   = (x**2 + y**2 + z**2)**0.5
    v   = (vx**2 + vy**2 + vz**2)**0.5

    eps = v**2 / 2 - mu / r # specific orbital energy [m^2/s^2]
    a   = -mu / 2 / eps     # semi-major axis [m]
    
    theta = np.arccos((x * vx + y * vy + z * vz) / r / v) # angle between r and v [rad]
    h     = r * v * np.sin(theta)                         # specific angular momentum [m^2/s]
    hx, hy, hz = y * vz - z * vy, z * vx - x * vz, x * vy - y * vx
    
    e  = (1 - h**2 / mu / a)**0.5                                      # eccentricity [rad]
    E  = np.arccos((1 - r / a) / e)                                    # eccentric anomaly [rad]
    M  = E - e * np.cos(E)                                             # mean anomaly [rad]
    i  = np.arccos(hz / h)                                             # inclination [rad]
    Ω  = (np.arctan2(hx, -hy) + 2 * np.pi) % (2 * np.pi)               # longitude of the ascending node [rad]
    nu = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))     # true anomaly [rad]
    ω  = np.arctan2(z / np.sin(i), x * np.cos(Ω) + y * np.sin(Ω)) - nu # argument of periapsis [rad]
    ϖ = (ω + Ω + 2 * np.pi) % (2 * np.pi)                              # longitude of periapsis [rad]
    L = (M + ϖ + 2 * np.pi) % (2 * np.pi)                              # mean longitude [rad]

    return e, a, i * 180 / np.pi, Ω * 180 / np.pi, ϖ * 180 / np.pi, L * 180 / np.pi
