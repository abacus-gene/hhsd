

class ParamTau():
    def __init__(self, tau_mean, tau_hpd_025, tau_hpd_975):
        self.mean = tau_mean
        self.hpd_025 = tau_hpd_025
        self.hpd_975 = tau_hpd_975
    
    def __str__(self):
        return f"Tau; mean: {self.mean}, hpd_025: {self.hpd_025}, hpd_975: {self.hpd_975}"

class ParamTheta():
    def __init__(self, theta_mean, theta_hpd_025, theta_hpd_975):
        self.mean = theta_mean
        self.hpd_025 = theta_hpd_025
        self.hpd_975 = theta_hpd_975

    def __str__(self):
        return f"Theta; mean: {self.mean}, hpd_025: {self.hpd_025}, hpd_975: {self.hpd_975}"

class ParamGDI():
    def __init__(self, gdi_mean, gdi_low, gdi_high):
        self.mean = gdi_mean
        self.low = gdi_low
        self.high = gdi_high

    def __str__(self):
        return f"GDI; mean: {self.mean}, lower: {self.low}, upper: {self.high}"