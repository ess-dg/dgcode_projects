from __future__ import print_function
import G4CustomPyGen
import Core.Units as Units
import Utils.NeutronMath
import math
from SansUtils.SansSphereGen import genSansSphereQR

class SansSphereGen(G4CustomPyGen.GenBase):
    """Generator using the SansSphereGen to generate neutrons following ideal SANS
       sphere-model distributions. Optionally, only neutrons in a given angular region
       can be selected."""

    def declare_parameters(self):
        self.addParameterDouble("zpos_mm",0.0)
        self.addParameterDouble("neutron_wavelength_aa",0.6)
        #self.addParameterDouble("sphere_radius_aa",60.0)
        self.addParameterDouble("sphere_radius_aa",10.0)
        self.addParameterDouble("polar_angle_min_degree",0.0)
        self.addParameterDouble("polar_angle_max_degree",180.0)
        self.addParameterDouble("fix_phi_degree",-1.0)

    def init_generator(self,gun):
        #NB: we generate q-values and relate to scattering angle theta by:
        #  q = (4pi/lambda) * sin(theta/2)

        wl = self.neutron_wavelength_aa*Units.angstrom
        sphere_radius = self.sphere_radius_aa*Units.angstrom
        self._k1 = self.neutron_wavelength_aa / ( self.sphere_radius_aa * 4.0 * math.pi )
        self._scatmin = max(0.0,self.polar_angle_min_degree*Units.degree)
        self._scatmax = min(math.pi,self.polar_angle_max_degree*Units.degree)
        assert self._scatmax>self._scatmin
        self._qrmin = math.sin( 0.5 * self._scatmin) / self._k1

        self._fixphi = self.fix_phi_degree * Units.degree if self.fix_phi_degree >= 0.0 else None
        gun.set_type('neutron')
        gun.set_position(0,0,self.zpos_mm*Units.mm)
        gun.set_energy(Utils.NeutronMath.neutronWavelengthToEKin(wl))
        self._nwarn = 10

    def generate_event(self,gun):
        while True:
            qr = genSansSphereQR(self.rand,self._qrmin)
            asinarg = self._k1 * qr
            if asinarg>1.0:
                #model breakdown...
                if self._nwarn:
                    print("SansSphereGen WARNING: Ignored out-of-range q value", end='')
                    self._nwarn -= 1
                    print(('' if self._nwarn else ' (suppressing further warnings)'))
                continue
            scat_angle = 2.0*math.asin(asinarg)
            if self._scatmin <= scat_angle <= self._scatmax:
                break
        s = math.sin(scat_angle)
        c = math.cos(scat_angle)
        if self._fixphi is None:
            phi = self.rand() * 2.0 * math.pi
        else:
            phi = self._fixphi
        gun.set_direction(s*math.cos(phi),s*math.sin(phi),c)
