"""Xylose -> I -> furfural, with F loss and bimolecular I+F loss.

Dimensionless time u=k_X*t. Parameters from Wang SI Table S9 at 160 C,
model [H3O+] = 0.10 M: alpha=k_I/k_X; beta=k_F/k_X;
gamma=k_FI*C_X0/k_X. Concentration states X,I,F normalized to inlet X.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, brentq

@dataclass(frozen=True)
class XyloseModel:
    alpha: float = 30.2348
    beta: float = 0.140586
    gamma_per_molar: float = 210.419
    k_x_per_second: float = 0.000259336

    def rates(self, y, concentration):
        x,i,f = y
        g = self.gamma_per_molar*concentration
        return np.array([-x, x-self.alpha*i-g*i*f,
                         self.alpha*i-self.beta*f-g*i*f])

    def jacobian(self,y,concentration):
        x,i,f=y; g=self.gamma_per_molar*concentration
        return np.array([[-1,0,0],[1,-self.alpha-g*f,-g*i],
                         [0,self.alpha-g*f,-self.beta-g*i]])

    def pfr(self,da,concentration):
        if da<0 or concentration<=0: raise ValueError('Require Da>=0 and C_X0>0')
        sol=solve_ivp(lambda t,y:self.rates(y,concentration),(0,da),
                      [1.,0.,0.],rtol=2e-10,atol=2e-12)
        if not sol.success: raise RuntimeError(sol.message)
        return sol.y[:,-1]

    def cstr(self,da,concentration,n):
        if int(n)!=n or n<1 or da<0 or concentration<=0: raise ValueError('Invalid arguments')
        h=da/n; g=self.gamma_per_molar*concentration
        x,i,f=1.,0.,0.
        for _ in range(n):
            xn=x/(1+h)
            if h==0: continue
            def i_of_fn(fn): return (i+h*xn)/(1+h*self.alpha+h*g*fn)
            def fun(fn): return fn*(1+h*self.beta+h*g*i_of_fn(fn))-f-h*self.alpha*i_of_fn(fn)
            fn=brentq(fun,0.,f+h*self.alpha*(i+h*xn)+1e-12,xtol=1e-14)
            i=i_of_fn(fn); x=xn; f=fn
        return np.array([x,i,f])

    def optimum(self,concentration,n=None,da_bounds=(0.05,20.)):
        fn=(lambda da: self.pfr(da,concentration)[2]) if n is None else (lambda da: self.cstr(da,concentration,n)[2])
        r=minimize_scalar(lambda da:-fn(da),bounds=da_bounds,method='bounded',options={'xatol':2e-9})
        if not r.success: raise RuntimeError('Optimization failed')
        return float(r.x),float(-r.fun)

    def phi(self,concentration):
        """Adjoint-equivalent first-order backward-Euler perturbation at PFR optimum.
        e'(u)=J(y)e+0.5 Da J(y)f(y); e(0)=0; y_N=y_PFR+e(Da)/N+O(N^-2).
        For optimized F yield: Phi=-e_F / F_PFR.
        """
        da,yopt=self.optimum(concentration)
        def rhs(u,v):
            y=v[:3]; e=v[3:]; ff=self.rates(y,concentration); jj=self.jacobian(y,concentration)
            return np.r_[ff,jj@e+0.5*da*(jj@ff)]
        sol=solve_ivp(rhs,(0,da),[1.,0.,0.,0.,0.,0.],rtol=2e-10,atol=2e-12)
        if not sol.success: raise RuntimeError(sol.message)
        return float(-sol.y[5,-1]/yopt)

    def stage_requirement(self,concentration,retention=0.95,max_stages=200):
        if not 0<retention<=1: raise ValueError('retention must lie in (0,1]')
        _,ref=self.optimum(concentration)
        for n in range(1,max_stages+1):
            if self.optimum(concentration,n)[1]>=retention*ref:
                return n
        return None

    def phi_sensitivity(self,concentration,parameter='gamma_per_molar',rel_step=0.002):
        """Central derivative of Phi relative to log(parameter), dimensionless."""
        from dataclasses import replace
        if parameter not in ('alpha','beta','gamma_per_molar'):
            raise ValueError('Choose alpha, beta, or gamma_per_molar')
        v=getattr(self,parameter)
        plus=replace(self,**{parameter:v*np.exp(rel_step)})
        minus=replace(self,**{parameter:v*np.exp(-rel_step)})
        return (plus.phi(concentration)-minus.phi(concentration))/(2*rel_step)
