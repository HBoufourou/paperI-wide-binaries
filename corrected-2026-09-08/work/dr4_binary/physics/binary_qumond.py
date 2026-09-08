"""Nonlinear QUMOND forces on TWO prescribed Plummer sources.

Dimensionless G=Mtotal=a0=1. Both mean forces are evaluated independently.
This module does not construct orbits, equilibrium populations, or Gaia tests.
"""
from functools import lru_cache
import time
import numpy as np

LEVELS = {
    "coarse": dict(log_width=.6, radial_order=6, nmu=24, nphi=48),
    "fine": dict(log_width=.45, radial_order=8, nmu=40, nphi=80),
    "reference": dict(log_width=.3, radial_order=10, nmu=64, nphi=128),
}


def nu_parts(y):
    """nu-1 and dnu/dy; rationalization preserves the high-field correction."""
    y = np.asarray(y, float)
    if np.any(y <= 0):
        raise ValueError("quadrature hit a zero of the Newtonian field")
    root = np.sqrt(.25+1/y)
    return 1/(y*(root+.5)), -.5/(y*y*root)


@lru_cache(maxsize=12)
def gauss(n):
    return np.polynomial.legendre.leggauss(n)


def logarithmic_rule(rmin, rmax, landmarks, width, order):
    boundaries = sorted(set([np.log(rmin), np.log(rmax)] +
                            [np.log(x) for x in landmarks if rmin < x < rmax]))
    node, weight = gauss(order)
    xs, ws = [], []
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        edges = np.linspace(a, b, int(np.ceil((b-a)/width))+1)
        for left, right in zip(edges[:-1], edges[1:]):
            xs.append((left+right)/2+(right-left)/2*node)
            ws.append((right-left)/2*weight)
    return np.exp(np.concatenate(xs)), np.concatenate(ws)


@lru_cache(maxsize=12)
def angular_rule(nmu, nphi):
    mu, w = gauss(nmu)
    phi = np.arange(nphi)*2*np.pi/nphi
    mm, pp = np.meshgrid(mu, phi, indexing="ij")
    d = np.stack((np.sqrt(1-mm*mm)*np.cos(pp), np.sqrt(1-mm*mm)*np.sin(pp), mm), axis=-1)
    weights = np.broadcast_to(w[:, None]*2*np.pi/nphi, mm.shape)
    return d.reshape(-1, 3), weights.ravel()


def binary_configuration(q=1., separation=1., theta_degrees=45., epsilon=.00125, external=1.):
    if q <= 0 or separation <= 0 or epsilon <= 0 or external < 0:
        raise ValueError("positive physical scales required")
    masses = np.array([1/(1+q), q/(1+q)])
    theta = np.deg2rad(theta_degrees)
    d = separation*np.array([np.sin(theta), 0., np.cos(theta)])
    positions = np.array([-masses[1]*d, masses[0]*d])
    return masses, positions, epsilon*np.sqrt(masses), np.array([0., 0., external])


def newton_plummer_force_magnitude(masses, b, separation, swapped=False):
    """Derivative of the exact angularly averaged two-Plummer interaction.

    This 1D integral is distinct from the main multicentre 3D quadrature.
    No effective softening sqrt(b1^2+b2^2) is introduced.
    """
    if swapped:
        b = b[::-1]
    d = float(separation)
    landmarks = [b[0]/2, b[0], 2*b[0], d/2, d, 2*d]
    for offset in [b[1], 10*b[1], 100*b[1]]:
        landmarks.extend([d-offset, d+offset])
    r, wt = logarithmic_rule(min(b)*1e-12, max(d, max(b))*1e6,
                             landmarks, .12, 16)
    plus = np.sqrt((r+d)**2+b[1]**2)
    minus = np.sqrt((r-d)**2+b[1]**2)
    derivative = 2*((d+r)/plus+(d-r)/minus)/(plus+minus)**2
    shell = 3*b[0]**2*r**3/(r*r+b[0]**2)**2.5
    return float(np.prod(masses)*np.sum(wt*shell*derivative))


def _source_geometry(x, mass, centre, b):
    r = x-centre
    r2 = np.sum(r*r, axis=1)
    k = mass/(r2+b*b)**1.5
    h = 3*mass/(r2+b*b)**2.5
    u = k[:, None]*r
    lap = h*b*b
    return r, r2, k, h, u, lap


def _h_times(geometry, vector):
    r, r2, k, h, u, lap = geometry
    return k[:, None]*vector-(h*np.sum(r*vector, axis=1))[:, None]*r


def _phantom_source(u, geometries):
    eta, dn = nu_parts(np.linalg.norm(u, axis=1))
    lap = sum(g[5] for g in geometries)
    hu = sum(_h_times(g, u) for g in geometries)
    s = eta*lap+dn/np.linalg.norm(u, axis=1)*np.sum(u*hu, axis=1)
    return s, eta


def mean_forces(masses, positions, b, external_vector, *, level="fine",
                outer_factor=128., inner_factor=1e-12, include_newton=False):
    """Return independent mean forces, with gradient-energy/source diagnostics.

    ``positions`` and ``external_vector`` permit general common orientations.
    The supplied source masses must sum to one. No CM projection is applied.
    """
    start = time.perf_counter()
    masses, positions, b, ext = [np.asarray(x, float) for x in
                                (masses, positions, b, external_vector)]
    if masses.shape != (2,) or positions.shape != (2, 3) or b.shape != (2,) or ext.shape != (3,):
        raise ValueError("invalid two-source input shapes")
    if np.any(masses <= 0) or np.any(b <= 0) or abs(masses.sum()-1) > 1e-12:
        raise ValueError("positive masses/softenings and Mtotal=1 required")
    settings = LEVELS[level].copy() if isinstance(level, str) else dict(level)
    dvec = positions[1]-positions[0]
    separation = np.linalg.norm(dvec)
    if separation == 0:
        raise ValueError("nonzero separation required")
    e3 = dvec/separation
    side = ext-(ext@e3)*e3
    if np.linalg.norm(side) < 1e-10*max(1., np.linalg.norm(ext)):
        coordinate = np.eye(3)[np.argmin(np.abs(e3))]
        side = coordinate-(coordinate@e3)*e3
    e1 = side/np.linalg.norm(side)
    basis = np.column_stack((e1, np.cross(e3, e1), e3))
    directions, womega = angular_rule(settings["nmu"], settings["nphi"])
    directions = directions@basis.T
    anomaly_energy = np.zeros((2, 3))
    anomaly_source = np.zeros((2, 3))
    raw_newton_energy = np.zeros((2, 3))
    raw_newton_source = np.zeros((2, 3))
    point_count = 0
    for centre_index in [0, 1]:
        rmin = inner_factor*np.sqrt(masses[centre_index])
        rmax = outer_factor*max(1., separation)
        scale = b[centre_index]
        radii, wr = logarithmic_rule(rmin, rmax,
                    [scale/2, scale, 2*scale, separation/2, separation, 2*separation, 1.],
                    settings["log_width"], settings["radial_order"])
        for first in range(0, len(radii), 8):
            radial = radii[first:first+8]
            x = positions[centre_index]+(radial[:, None, None]*directions[None, :, :]).reshape(-1, 3)
            geometries = [_source_geometry(x, masses[i], positions[i], b[i]) for i in [0,1]]
            radii2 = [g[1] for g in geometries]
            opposite = 1-centre_index
            partition = radii2[opposite]**3/(radii2[0]**3+radii2[1]**3)
            volume = (wr[first:first+8, None]*radial[:, None]**3*womega[None, :]).ravel()
            weight = volume*partition/(4*np.pi)
            utotal = geometries[0][4]+geometries[1][4]-ext
            stotal, etatotal = _phantom_source(utotal, geometries)
            for i in [0, 1]:
                ui_external = geometries[i][4]-ext
                ssingle, etasingle = _phantom_source(ui_external, [geometries[i]])
                da = etatotal[:, None]*utotal-etasingle[:, None]*ui_external
                anomaly_energy[i] -= np.sum(weight[:, None]*_h_times(geometries[i], da), axis=0)
                anomaly_source[i] += np.sum((weight*(stotal-ssingle))[:, None]*geometries[i][4], axis=0)
                if include_newton:
                    raw_newton_energy[i] -= np.sum(weight[:, None]*_h_times(geometries[i], geometries[1-i][4]), axis=0)
                    # The own Newtonian self-force is analytically zero.
                    raw_newton_source[i] += np.sum((weight*geometries[1-i][5])[:, None]*geometries[i][4], axis=0)
            point_count += len(x)
    newton_magnitude = newton_plummer_force_magnitude(masses, b, separation)
    newton = np.array([newton_magnitude*e3, -newton_magnitude*e3])
    force = newton+anomaly_energy
    source_force = newton+anomaly_source
    force_scale = np.mean(np.linalg.norm(force, axis=1))
    result = dict(masses=masses.tolist(), positions=positions.tolist(), b=b.tolist(),
                  external_newtonian=ext.tolist(), separation=float(separation), settings=settings,
                  outer_factor=float(outer_factor), inner_factor=float(inner_factor),
                  force=force.tolist(), acceleration=(force/masses[:, None]).tolist(),
                  relative_acceleration=(force[1]/masses[1]-force[0]/masses[0]).tolist(),
                  source_force=source_force.tolist(), newton_force=newton.tolist(),
                  cm_relative_residual=float(np.linalg.norm(force.sum(axis=0))/force_scale),
                  source_energy_relative=float(np.max(np.linalg.norm(force-source_force, axis=1)/np.linalg.norm(force, axis=1))),
                  quadrature_points=point_count, elapsed_seconds=time.perf_counter()-start)
    if include_newton:
        result.update(raw_newton_energy=raw_newton_energy.tolist(), raw_newton_source=raw_newton_source.tolist(),
                      newton_swapped_magnitude=newton_plummer_force_magnitude(masses,b,separation,swapped=True))
    return result


def solve_case(q=1., separation=1., theta_degrees=45., epsilon=.00125,
               external=1., level="fine", **kwargs):
    result = mean_forces(*binary_configuration(q,separation,theta_degrees,epsilon,external),
                         level=level, **kwargs)
    result["case"] = dict(q=q, separation=separation, theta_degrees=theta_degrees,
                          epsilon=epsilon, external=external, level=level)
    return result


if __name__ == "__main__":
    import json
    print(json.dumps(solve_case(include_newton=True), indent=2), flush=True)
