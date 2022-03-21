from mrchem_input_generator import MRChemInputGenerator, Molecule, RecursiveNamespace


class EnergyCalculation(MRChemInputGenerator):
    """Simple class for auto-generating an MRChem energy calculation.
    Options can be overwritten by using the interface to MRChemInputGenerator."""
    def __init__(self, molecule=None, fname='energy.inp',
                 world_prec=1.0e-4, method='lda', kain_scf=6,
                 guess_type='sad_dz', localize=True, **kwargs):
        super().__init__(**kwargs)
        self.fname = fname

        # Add some sections
        self.add_input_section('Molecule', 'WaveFunction', 'SCF')

        # Load Molecule
        if molecule is not None:
            assert isinstance(molecule, Molecule), 'Error reading molecule: You need to pass a Molecule object.'
            self.molecule = molecule
            self.input.Molecule.coords = self.molecule.to_string_format()
            self.input.Molecule.charge = self.molecule.charge
            self.input.Molecule.multiplicity = self.molecule.multiplicity

        self.input.world_prec = world_prec
        self.input.SCF.kain = kain_scf
        self.input.SCF.guess_type = guess_type
        self.input.SCF.localize = localize
        self.input.WaveFunction.method = method

        self.defaults = self.get_defaults()

    def get_defaults(self):
        """Collect default values for the current input file."""
        i = self.to_dict(self.input)
        return RecursiveNamespace(**{
            key: val for key, val in self.reference.defaults.items() if key in i.keys()
        })

    def write(self):
        """Write input to file."""
        with open(self.fname, 'w') as f:
            f.write(str(self))


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

        self.defaults = self.get_defaults()


class MagneticResponseCalculation(EnergyCalculation):
    """Simple class for auto-generating an MRChem magnetic response calculation.
    Inherits from EnergyCalculation.
    Options can be overwritten by using the interface to MRChemInputGenerator."""
    def __init__(self, kain_rsp=6, magnetizability=True, nmr_shielding=True,
                 nuclei=None, nuclear_specific=False, **kwargs):
        super().__init__(**kwargs)
        self.fname = 'magnetic_rsp.inp'

        if nuclei is None and self.molecule is not None:
            nuclei = [i for i in range(self.molecule.n_atoms)]

        self.add_input_section('Properties', 'Response', 'NMRShielding')
        self.input.Properties.magnetizability = magnetizability
        self.input.Properties.nmr_shielding = nmr_shielding
        self.input.NMRShielding.nuclear_specific = nuclear_specific
        self.input.NMRShielding.nucleus_k = nuclei
        self.input.Response.kain = kain_rsp

        self.defaults = self.get_defaults()


class EnergyZORACalculation(EnergyCalculation):
    """Simple class for auto-generating an MRChem energy calculation with ZORA activated.
    Options can be overwritten by using the interface to MRChemInputGenerator."""
    def __init__(self, light_speed=-1.0, include_nuclear=True, include_coulomb=True, include_xc=True, **kwargs):
        super().__init__(**kwargs)
        self.fname = 'energy_zora.inp'

        self.add_input_section('ZORA')
        self.input.WaveFunction.relativity = True
        self.input.ZORA.light_speed = light_speed
        self.input.ZORA.include_nuclear = include_nuclear
        self.input.ZORA.include_coulomb = include_coulomb
        self.input.ZORA.include_xc = include_xc

        self.defaults = self.get_defaults()


if __name__ == '__main__':
    mol = Molecule.from_xyzfile('NH3O.xyz')
    e = EnergyCalculation(molecule=mol)
    rsp_e = ElectricResponseCalculation(molecule=mol)
    rsp_m = MagneticResponseCalculation(molecule=mol)
    e_zora = EnergyZORACalculation(molecule=mol)

    e_zora.input.ZORA.include_xc = False
    e.write()
    e_zora.write()
    rsp_e.write()
    rsp_m.write()
