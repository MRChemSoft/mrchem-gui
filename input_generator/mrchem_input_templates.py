from mrchem_input_generator import MRChemInputGenerator, PrettyDict


class EnergyCalculation(MRChemInputGenerator):
    """Simple class for auto-generating an MRChem energy calculation.
    Options can be overwritten by using the interface to MRChemInputGenerator."""
    def __init__(self, xyzfile=None, fname='energy.inp', world_prec=1e-4, method='lda', kain_scf=6, charge=0,
                 mult=1, unit='angstrom', guess_type='sad_dz', localize=True):
        super().__init__(hide_defaults=True)
        self.fname = fname
        self.xyzfile = xyzfile

        self.add_input_section('Molecule', 'WaveFunction', 'SCF')

        self.input.world_prec = world_prec
        self.input.world_unit = unit
        self.input.Molecule.charge = charge
        self.input.Molecule.multiplicity = mult

        if self.xyzfile is not None:
            self.xyzfile = xyzfile
            self._load_xyzfile()
            self.input.Molecule.coords = self.coords
        self.input.SCF.kain = kain_scf
        self.input.SCF.guess_type = guess_type
        self.input.SCF.localize = localize
        self.input.WaveFunction.method = method

        self._set_defaults()

    def _load_xyzfile(self):
        """Load XYZ file to correct format for an MRChem JSON input file."""
        with open(self.xyzfile) as f:
            lines = [line.strip() for line in f.readlines()]
        self.n_atoms = int(lines.pop(0))
        self.xyz_comment = lines.pop(0)
        self.coords = '\n'.join(lines)

    def _set_defaults(self):
        """Collect default values for the current input file."""
        self.defaults = PrettyDict(**{
            section: val for section, val in self.to_dict(self._defaults).items() if section in self.to_dict(self.input).keys()
        })

    def write(self):
        """Write input to file."""
        with open(self.fname, 'w') as f:
            f.write(str(self.generate_without_defaults()))


class ElectricResponseCalculation(EnergyCalculation):
    """Simple class for auto-generating an MRChem electric response calculation.
    Inherits from EnergyCalculation.
    Options can be overwritten by using the interface to MRChemInputGenerator."""
    def __init__(self, kain_rsp=6, polarizability=True, quadrupole_moment=True, frequencies=None,
                 field_strength=None, **kwargs):
        super().__init__(**kwargs)
        self.fname = 'electric_rsp.inp'

        if field_strength is None:
            field_strength = [0.001, 0.001, 0.001]
        if frequencies is None:
            frequencies = [0.0]

        self.add_input_section('ExternalFields', 'Response', 'Polarizability', 'Properties')

        self.input.Properties.quadrupole_moment = quadrupole_moment
        self.input.Properties.polarizability = polarizability
        self.input.Polarizability.frequency = frequencies
        self.input.ExternalFields.electric_field = field_strength
        self.input.Response.kain = kain_rsp

        self._set_defaults()


class MagneticResponseCalculation(EnergyCalculation):
    """Simple class for auto-generating an MRChem magnetic response calculation.
    Inherits from EnergyCalculation.
    Options can be overwritten by using the interface to MRChemInputGenerator."""
    def __init__(self, kain_rsp=6, magnetizability=True, nmr_shielding=True,
                 nuclei=None, nuclear_specific=False, **kwargs):
        super().__init__(**kwargs)
        self.fname = 'magnetic_rsp.inp'

        if nuclei is None and self.xyzfile is not None:
            nuclei = [i for i, atom in enumerate(self.coords.split('\n'))]

        self.add_input_section('Properties', 'Response', 'NMRShielding')
        self.input.Properties.magnetizability = magnetizability
        self.input.Properties.nmr_shielding = nmr_shielding
        self.input.NMRShielding.nuclear_specific = nuclear_specific
        self.input.NMRShielding.nucleus_k = nuclei
        self.input.Response.kain = kain_rsp

        self._set_defaults()


if __name__ == '__main__':
    print('Generating sample input files (energy, electric response, and magnetic response')
    e = EnergyCalculation(xyzfile='NH3O.xyz')
    rsp_e = ElectricResponseCalculation(xyzfile='NH3O.xyz')
    rsp_m = MagneticResponseCalculation(xyzfile='NH3O.xyz')

    e.write()
    rsp_e.write()
    rsp_m.write()
