import numpy as np

Msun = 1.9884e30 # [kg]
G    = 6.6743e-11 # [m^3/kg/s^2]

def kep2cart(e, a, i, Ω, m, ϖ=None, ω=None, L=None, M=None):
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
        Total mass (central body + orbiting object) [kg]
    ϖ: float
        Longitude of periapsis [deg]; used if ω is None
    ω: float
        Argument of periapsis [deg]
    L: float
        Mean longitude [deg] at given epoch; used if M is None
    M: float
        Mean anomaly [deg] at given epoch

    Returns
    -------
    (x, y, z), (vx, vy, vz): tuples
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
        ω = (ϖ - Ω) * np.pi / 180 # argument of periapsis [rad], 0 < ω < 2*pi
    else: # provided ω
        ω = ω * np.pi / 180
    
    i = i * np.pi / 180
    Ω = Ω * np.pi / 180

    if M is None:
        M = L - ϖ # mean anomaly [deg], 0 < M < 360

    E  = ecc_anomaly(e, M) * np.pi / 180 # eccentric anomaly [rad], 0 < E < 2*pi
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), \
                        np.sqrt(1 - e) * np.cos(E / 2)) # true anomaly [rad], 0 < nu < 2*pi

    # heliocentric coordinates [m] in orbital plane, z0 = 0
    r0 = a * (1 - e * np.cos(E))
    x0 = r0 * np.cos(nu) # = a * (np.cos(E) - e)
    y0 = r0 * np.sin(nu) # = a * (1 - e**2)**0.5 * np.sin(E)

    # ecliptic coordinates [m]
    x =   ( np.cos(ω) * np.cos(Ω) - np.sin(ω) * np.sin(Ω) * np.cos(i)) * x0 \
        + (-np.sin(ω) * np.cos(Ω) - np.cos(ω) * np.sin(Ω) * np.cos(i)) * y0
    y =   ( np.cos(ω) * np.sin(Ω) + np.sin(ω) * np.cos(Ω) * np.cos(i)) * x0 \
        + (-np.sin(ω) * np.sin(Ω) + np.cos(ω) * np.cos(Ω) * np.cos(i)) * y0
    z =   np.sin(ω) * np.sin(i) * x0 + np.cos(ω) * np.sin(i) * y0

    # find velocity in ecliptic coordinates [m/s]
    mu = G * m # [m^3/s^2]
    v  = np.sqrt(mu / a * (1 + e * np.cos(E)) / (1 - e * np.cos(E))) # orbital speed [m/s]
    # INCORRECT
    #θ = np.arcsin(np.sqrt((1 - e**2) / (1 - e**2 * np.cos(E)**2))) # angle between r and v [rad]
    #vx0 = v * np.cos(nu + θ) # velocity in orbital plane [m/s]
    #vy0 = v * np.sin(nu + θ)
    # CORRECT
    vx0 = -v / np.sqrt(1 - e**2 * np.cos(E)**2) * np.sin(E) # velocity in orbital plane [m/s]
    vy0 =  v / np.sqrt(1 - e**2 * np.cos(E)**2) * np.sqrt(1 - e**2) * np.cos(E) 

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

def cart2kep(x, y, z, vx, vy, vz, m):
    """
    Convert cartesian state vectors (r [m] and v [m/s]) of a given epoch 
    in the ecliptic coordinates to Keplerian orbital elements. The x-axis 
    points toward the vernal equinox (in the solar system case).

    Parameters
    ----------
    x, y, z: float
        Position [m] measured in the initial rest frame of the central body
    vx, vy, vz: float
        Velocity [m/s] measured in the initial rest frame of the central body
    m: float
        Total mass (central body + orbiting object) [kg]

    Returns
    -------
    e, a, i, Ω, ϖ, L: tuple
        Eccentricity (0 < e < 1), semi-major axis a [m], inclination i [deg] 
        (0 < i < 180), longitude of the ascending node Ω [deg] (0 < Ω < 360), 
        longitude of periapsis ϖ [deg] (0 < ϖ < 360), and mean longitude L [deg]
        at given epoch (0 < L < 360)

    Notes
    -----
    (1) Since the inverse trig functions have limited ranges (0 < np.arccos < pi, 
        -pi/2 < np.arcsin < pi/2, -pi < np.arctan2 < pi), one must be careful in
        the calculation of the eccentric anomaly and the angle between the position
        and velocity state vectors to account for all possible values
    (2) Due to (1), we cannot use E  = np.arccos((1 - r / a) / e) to caluclate the
        eccentric anomaly because it only has a range of 0 < E < pi.
    (3) cosθ, the cosine between the state vectors (x, y, z) and (vx, vy, vz), 
        should be positive when pi/2 < E < pi (quadrant II) or 3*pi/2 < E < 2*pi 
        (quadrant IV), and negative otherwise.
    """
    mu   = G * m # [m^3/s^2]
    r    = (x**2 + y**2 + z**2)**0.5
    v    = (vx**2 + vy**2 + vz**2)**0.5
    cosθ = (x * vx + y * vy + z * vz) / r / v # cosine between (x, y, z) and (vx, vy, vz)
    
    eps  = v**2 / 2 - mu / r # specific orbital energy [m^2/s^2]
    a    = -mu / 2 / eps     # semi-major axis [m]
    
    # h is conserved, so sinθ should stay positive
    h    = r * v * np.sqrt(1 - cosθ**2) # specific angular momentum [m^2/s]
    hx, hy, hz = y * vz - z * vy, z * vx - x * vz, x * vy - y * vx
    
    e  = (1 - h**2 / mu / a)**0.5                        # eccentricity [rad]
    i  = np.arccos(hz / h)                               # inclination [rad], 0 < i < pi
    Ω  = (np.arctan2(hx, -hy) + 2 * np.pi) % (2 * np.pi) # longitude of the ascending node [rad], 0 < Ω < 2*pi
    
    p  = a * (1 - e**2) # semi-latus rectum [m]
    # INCORRECT - LIMITED ANGLE RANGE
    #E  = np.arccos((1 - r / a) / e)                                                     # eccentric anomaly [rad], 0 < E < pi
    #nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2)) # true anomaly [rad], 0 < nu < pi
    # CORRECT
    nu = (np.arctan2(np.sqrt(p / mu) * r * v * cosθ, p - r) + 2 * np.pi) % (2 * np.pi)    # true anomaly [rad], 0 < nu < 2*pi
    E  = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2)) # eccentric anomaly [rad], 0 < E < 2*pi
    M  = E - e * np.sin(E)                                                                # mean anomaly [rad], 0 < M < 2*pi

    ω  = np.arctan2(z / np.sin(i), x * np.cos(Ω) + y * np.sin(Ω)) - nu # argument of periapsis [rad]
    ϖ = (ω + Ω + 2 * np.pi) % (2 * np.pi)                              # longitude of periapsis [rad]
    L = (M + ϖ + 2 * np.pi) % (2 * np.pi)                              # mean longitude [rad]

    return e, a, i * 180 / np.pi, Ω * 180 / np.pi, ϖ * 180 / np.pi, L * 180 / np.pi

def result2kep(input_file, name, central_body, snapshot=-1):
    """
    Calculate orbital elements from simuation results

    Parameters
    ----------
    input_file: str
        Path to simulation result
    name: str
        Name of the orbiting body whose orbital elements are to be calculated
    central_body: str
        Name of the central body
    snapshot: int, optional
        Index of the snapshot used; defaut: -1 (last)

    Returns
    -------
    e, a, i, Ω, ϖ, L, m: tuple
    """
    result = np.load(input_file, allow_pickle=True).item()

    i = result['id_index'][name]
    X, V = result['x'][snapshot, i], result['v'][snapshot, i]
    m = result['kepler'][name][6]
    m0 = result['kepler'][central_body][6]

    return cart2kep(*X, *V, m + m0)

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        assert len(sys.argv) >= 4
        print('(e, a, i, Ω, ϖ, L) = ', result2kep(*sys.argv[1:4]))
