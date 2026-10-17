"""Exploratory local-circular-arc model of a two-front liquid plug.

Millimetres, radians, areas in mm², volumes in uL (1 mm³ = 1 uL).
This is NOT calibrated against measured wall profiles or trajectories.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline
from scipy.linalg import null_space
from scipy.optimize import brentq, root, minimize


class PlugModel:
    def __init__(self, theta_deg=40, width_mm=1., amplitude_mm=.05,
                 period_mm=.30, depth_mm=.38):
        self.th = np.deg2rad(theta_deg)
        self.w0 = width_mm
        self.a = amplitude_mm
        self.k = 2 * np.pi / period_mm
        self.depth = depth_mm
        x = np.linspace(-4, 4, 160001)
        q = np.sqrt(1 + (self.dw(x) / 2)**2)
        self.S = CubicSpline(x, cumulative_trapezoid(q, x, initial=0))

    def w(self, x):
        return self.w0 + self.a * np.cos(self.k * np.asarray(x))

    def dw(self, x):
        return -self.a * self.k * np.sin(self.k * np.asarray(x))

    def ddw(self, x):
        return -self.a * self.k**2 * np.cos(self.k * np.asarray(x))

    def B(self, x):
        return self.w0 * np.asarray(x) + self.a / self.k * np.sin(self.k * np.asarray(x))

    def alpha(self, x):
        return np.arctan(self.dw(x) / 2)

    def theta(self, x, side):
        return self.th + (1 if side == "right" else -1) * self.alpha(x)

    @staticmethod
    def f(t):
        t = np.asarray(t)
        return np.where(abs(t - np.pi / 2) < 1e-5,
                        -(np.pi / 2 - t) / 6,
                        (2*t - np.pi + np.sin(2*t)) / (8*np.cos(t)**2))

    @staticmethod
    def fp(t):
        return .5 + 2 * PlugModel.f(t) * np.tan(t)

    @staticmethod
    def g(t):
        return (np.pi - 2*t) / (2 * np.cos(t))

    def volume(self, left, right):
        return self.depth * (self.B(right) - self.B(left) +
            self.w(left)**2 * self.f(self.theta(left, "left")) +
            self.w(right)**2 * self.f(self.theta(right, "right")))

    def p(self, x, side):
        return 2 * np.cos(self.theta(x, side)) / self.w(x)

    def pprime(self, x, side):
        w, dw, dd = self.w(x), self.dw(x), self.ddw(x)
        t = self.theta(x, side)
        sign = 1 if side == "right" else -1
        tp = sign * (dd / 2) / (1 + (dw / 2)**2)
        return 2 * (-np.sin(t)*tp*w - np.cos(t)*dw) / w**2

    def symmetric_state(self, s):
        return np.array([-s, s,
                         float(self.theta(-s, "left")),
                         float(self.theta(s, "right"))])

    def energy_area(self, z):
        left, right, tl, tr = z
        wl, wr = self.w(np.array([left, right]))
        area = (self.B(right)-self.B(left) + wl**2*self.f(tl) +
                wr**2*self.f(tr))
        energy = (wl*self.g(tl) + wr*self.g(tr) -
                  2*np.cos(self.th)*(self.S(right)-self.S(left)))
        return float(energy), float(area)

    def gradients(self, z):
        left, right, tl, tr = z
        x = np.array([left, right])
        t = np.array([tl, tr])
        w, dw = self.w(x), self.dw(x)
        q = np.sqrt(1+(dw/2)**2)
        fp = self.fp(t)
        gp = -2*np.cos(t)*fp
        grad_e = np.r_[dw*self.g(t)+np.array([1,-1])*2*np.cos(self.th)*q, w*gp]
        grad_a = np.r_[np.array([-1,1])*w+2*w*dw*self.f(t), w*w*fp]
        return grad_e, grad_a

    def constrained_stability(self, z):
        """Eigenvalues of constrained Hessian of E + lambda*A."""
        left, right, tl, tr = z
        x = np.array([left, right])
        t = np.array([tl, tr])
        w, dw, dd = self.w(x), self.dw(x), self.ddw(x)
        signs = np.array([-1.,1.])
        q = np.sqrt(1+(dw/2)**2)
        qp = dw*dd/(4*q)
        multiplier = 2*np.cos(tr)/w[1]
        f, fp = self.f(t), self.fp(t)
        fpp = 2*fp*np.tan(t) + 2*f/np.cos(t)**2
        g = self.g(t)
        gp = -2*np.cos(t)*fp
        gpp = 2*np.sin(t)*fp - 2*np.cos(t)*fpp
        grad_e, grad_a = self.gradients(z)
        H = np.zeros((4,4))
        for i in range(2):
            H[i,i] = (dd[i]*g[i] - signs[i]*2*np.cos(self.th)*qp[i] +
                      multiplier*(signs[i]*dw[i] +
                         2*(dw[i]**2+w[i]*dd[i])*f[i]))
            H[i,2+i] = H[2+i,i] = (dw[i]*gp[i] +
                multiplier*2*w[i]*dw[i]*fp[i])
            H[2+i,2+i] = w[i]*gpp[i]+multiplier*w[i]**2*fpp[i]
        basis = null_space(grad_a[None,:])
        eigenvalues = np.linalg.eigvalsh(basis.T @ H @ basis)
        residual = np.linalg.norm(grad_e+multiplier*grad_a)
        return eigenvalues, residual

    def critical_and_post(self):
        """A particular symmetric branch and one positive-direction landing."""
        s = brentq(lambda v: float(self.pprime(v, "right")), .67, .735)
        sf = brentq(lambda v: float(self._sym_area_grad_right(v)), .62, .67)
        initial = self.symmetric_state(s)
        target = self.volume(-s, s)
        def eq(pos):
            l, r = pos
            return [float(self.p(l, "left")-self.p(r, "right")),
                    float(self.volume(l, r)-target)]
        solution = root(eq, [-.44,.74], tol=1e-12)
        if max(abs(v) for v in eq(solution.x)) > 1e-8:
            raise RuntimeError("Post-jump pressure/volume root failed")
        l, r = solution.x
        final = np.array([l, r, float(self.theta(l, "left")),
                         float(self.theta(r, "right"))])
        eig, residual = self.constrained_stability(final)
        if min(eig) < 1e-8 or residual > 1e-7:
            raise RuntimeError("Landing state not independently stable")
        return dict(critical_volume_uL=float(target),
                    symmetric_fold_volume_uL=float(self.volume(-sf,sf)),
                    half_length_at_instability_mm=float(s),
                    before=initial, after=final,
                    predicted_length_drop_um=float((2*s-(r-l))*1000),
                    predicted_center_shift_um=float((l+r)*500),
                    min_landing_eigenvalue=float(min(eig)),
                    energy_drop=float(self.energy_area(initial)[0]-
                                      self.energy_area(final)[0]))

    def _sym_area_grad_right(self, s):
        # Right coordinate derivative with full circular-arc area correction.
        x = s
        w, dw, dd = self.w(x), self.dw(x), self.ddw(x)
        theta = self.theta(x, "right")
        dt = (dd/2)/(1+(dw/2)**2)
        return w + 2*w*dw*self.f(theta) + w*w*self.fp(theta)*dt


def independent_minimization(model, target_volume_uL, seed):
    """Full four-variable constrained optimization, not just a pressure root."""
    target_area = target_volume_uL/model.depth
    result = minimize(lambda z: model.energy_area(z)[0], seed,
        jac=lambda z: model.gradients(z)[0], method="SLSQP",
        constraints=[dict(type="eq", fun=lambda z:model.energy_area(z)[1]-target_area,
                          jac=lambda z:model.gradients(z)[1])],
        bounds=[(-1.8,-.015),(.015,1.8),(.01,2.9),(.01,2.9)],
        options=dict(ftol=1e-12,maxiter=400))
    if not result.success:
        raise RuntimeError(result.message)
    return result.x
