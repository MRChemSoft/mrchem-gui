import copy
import json
import yaml
import requests
import numpy as np
import mendeleev
from mendeleev import element
from mendeleev.econf import ElectronicConfiguration
from types import SimpleNamespace


def couldbefloat(num):
    """Helper function"""
    try:
        float(num)
        return True
    except ValueError:
        return False


class RecursiveNamespace(SimpleNamespace):
    """Recursively convert nested dict to SimpleNamespaces, to achieve
    a dotted access to all keys. Inspired from

    https://dev.to/taqkarim/extending-simplenamespace-for-nested-dictionaries-58e8"""
    @staticmethod
    def map_to_ns(entry):
        """Map [dict] to [SimpleNamespace]"""
        if isinstance(entry, dict):
            return RecursiveNamespace(**entry)

    @staticmethod
    def map_to_dict(entry):
        """Map [SimpleNamespaces] to [dict]"""
        if isinstance(entry, RecursiveNamespace):
            return entry.__dict__

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        for key, val in kwargs.items():
            if isinstance(val, dict):
                setattr(self, key, RecursiveNamespace(**val))
            elif isinstance(val, list):
                setattr(self, key, val)

    def __str__(self):
        return json.dumps(MRChemInputGenerator.to_dict(self), indent=MRChemInputGenerator.json_indentation)


class MRChemInputGenerator:
    """
    Class for generating MRChem input files. All input blocks and keywords are fetched from the
    input reference `template.yml` at the Github repo for release version 1.

    Input blocks are represented as `SimpleNamespace`s, to allow for dot access to specific keywords.

    Example usage:

    ```python
    from mrchem_input_generator import MRChemInputGenerator

    m = MRChemInputGenerator(hide_defaults=True)
    m.add_input_section('Molecule', 'SCF', 'WaveFunction', 'MPI')
    m.input.SCF.kain = 10
    m.input.SCF.max_iter = 100
    m.input.WaveFunction.method = 'lda'
    m.input.WaveFunction.restricted = False
    m.input.Molecule.charge = 1
    m.input.Molecule.multiplicity = 4
    m.input.MPI.numerically_exact = True
    print(m.generate_without_defaults())

    with open('mrchem.inp', 'w') as f:
        f.write(str(mrc.generate()))
    ```

    which generates the following input file

    ```
    {
        "world_prec": null,
        "world_size": -1,
        "world_unit": "bohr",
        "world_origin": [
            0.0,
            0.0,
            0.0
        ],
        "Molecule": {
            "charge": 1,
            "multiplicity": 4
        },
        "SCF": {
            "kain": 10
        },
        "WaveFunction": {
            "method": "lda",
            "restricted": false
        },
        "MPI": {
            "numerically_exact": true
        }
    }
    ```


    """
    # _remote_template = 'https://raw.githubusercontent.com/MRChemSoft/mrchem/release/1.0/python/template.yml'
    remote_template = 'https://raw.githubusercontent.com/Andersmb/mrchem/new-zora/python/template.yml'
    json_indentation = 2

    def __init__(self, hide_defaults=True):
        self.hide_defaults = hide_defaults
        self.template = yaml.safe_load(requests.get(self.remote_template).text)

        self.input = RecursiveNamespace(**{})
        self._defaults = RecursiveNamespace(**self._parse_tmplt())

        # Add the top level 'world' keywords (always present in the input file)
        for key in [key for key in self.to_dict(self._defaults).keys() if 'world' in key]:
            self.input.__setattr__(key, self._defaults.__getattribute__(key))

    def __str__(self):
        """Pretty formatted string representation"""
        if not self.hide_defaults:
            return json.dumps(self.to_dict(self.input), indent=self.json_indentation)
        else:
            i = self.to_dict(self.input)
            d = self.to_dict(self._defaults)

            keys_to_delete = []
            filtered = copy.deepcopy(i)
            for section in i.keys():
                if isinstance(i[section], dict):
                    for key, keyval in i[section].items():
                        ref = d[section][key]
                        val = i[section][key]
                        if val == ref:
                            keys_to_delete.append((section, key))

            for section, key in keys_to_delete:
                del filtered[section][key]
                if filtered[section] == {}:
                    del filtered[section]

            return json.dumps(filtered, indent=self.json_indentation)

    @classmethod
    def from_json(cls, file):
        """Instantiate MRChemInputGenerator object from existing input file in JSON format."""
        with open(file) as f:
            j = json.loads(f.read())
        inp = cls()
        for key, val in j.items():
            inp.input.__setattr__(key, val)
        return inp

    def from_text(self, file):
        """Instantiate MRChemGenerator object from existing input file in text format.
        Needs to parse the text file into correct JSON.
        """
        raise NotImplementedError

    def _parse_tmplt(self) :
        """
        Read the MRChem input template, and generate a nested dict of all keywords
        and their defaults (if present), which can easily be turned in to a recursive namespace
        (i.e. to get dotted access to all keys).
        """
        # The top level 'world' keywords
        top_level = {
            key['name']: key['default'] if 'default' in key.keys() else None for key in self.template['keywords']
        }

        # The remaining sections and their keywords
        sections = {
            section['name']: {
                key['name'].lower(): key['default'] if 'default' in key.keys() else None for key in section['keywords']
            }
            for section in self.template['sections']
        }

        return {**top_level, **sections}

    def add_input_section(self, *sections):
        """Add input section."""
        for section in sections:
            if not hasattr(self._defaults, section):
                raise InvalidInputSection(section)
            self.input.__setattr__(
                section, RecursiveNamespace(**self._defaults.__getattribute__(section).__dict__)
            )

    def remove_input_section(self, *sections):
        """Delete input section."""
        for section in sections:
            if not hasattr(self._defaults, section):
                raise InvalidInputSection(section)
            self.input.__delattr__(section)

    def set_keyword(self, section, key, value):
        """Set/update a (key, val) pair in the passed section"""
        # Determine whether the value is float, int, bool, or str
        # Argparse produces strings by default from nargs option

        # Validate section and key by comparing to template.yml
        if not hasattr(self._defaults, section):
            raise InvalidInputSection(section)
        elif not hasattr(getattr(self._defaults, section), key):
            raise InvalidInputKeyword(key)

        # Check if section exists. If not add it
        if not hasattr(self.input, section):
            self.add_input_section(section)

        # Convert from RecursiveNamespace to dict
        d = self.to_dict(self.input)

        # Make sure the correct type is used
        # Order of tests:
        # float, +int, -int, float, True, False, str
        # TODO: There must be an easier way to get correct types!
        # What about lists?
        if couldbefloat(value) and '.' in value:
            d[section][key] = float(value)
        elif value.isdigit():
            d[section][key] = int(value)
        elif value.startswith('-') and value[1:].isdigit():
            d[section][key] = int(value)
        elif couldbefloat(value):
            d[section][key] = float(value)
        elif value.lower() == 'false':
            d[section][key] = False
        elif value.lower() == 'true':
            d[section][key] = True
        else:
            d[section][key] = value

        # Convert back to RecursiveNamespace
        self.input = RecursiveNamespace(**d)

    @staticmethod
    def to_dict(ns):
        """Recursively convert SimpleNamespaces to dictionaries."""
        d = {}
        for key, val in ns.__dict__.items():
            if isinstance(val, SimpleNamespace):
                d[key] = val.__dict__
            else:
                d[key] = val
        return d


class Molecule:
    """Simple class for storing information about the atomic/molecular system to be studied."""
    def __init__(self, xyz_file=None, charge=0, multiplicity=1, precision=10):
        self.xyz_file = xyz_file                              # XYZ file storing coordinates
        self.charge = charge                                  # Total charge
        self.multiplicity = multiplicity                      # Total spin multiplicity
        self.precision = precision                            # Number of decimals in coordinate
        self.n_unpaired_electrons = self.multiplicity - 1     # Number of unpaired electrons

        # Load XYZ file if given
        if self.xyz_file is not None:
            self.n_atoms, self.comment, self.symbols, self.coords = self.load_xyzfile(self.xyz_file)

    def __str__(self):
        """Pretty print coordinate section of XYZ file."""
        symbols = self.symbols.tolist()
        coords = self.coords.tolist()
        new = []
        for sym, c in zip(symbols, coords):
            x, y, z = c
            prettify = lambda c: f'{c:.{self.precision}f}' if c < 0 else f' {c:.{self.precision}f}'
            x = prettify(x)
            y = prettify(y)
            z = prettify(z)

            new.append(f'{sym[0]} {x} {y} {z}')
        return '\n'.join(new)

    def __eq__(self, other):
        """Check if two Molecules are identical (parameters related to calculation are checked)."""
        if not isinstance(other, Molecule):
            return False
        tests = [
            (self.coords == other.coords).all(),
            self.charge == other.charge,
            self.multiplicity == other.multiplicity
        ]

        return all(tests)

    def to_xyz_format(self):
        """Return coordinates in XYZ file format."""
        return f'{self.n_atoms}\n' + f'{self.comment}\n' + str(self)

    def to_string_format(self, delimiter='\n'):
        """Return coordinates in <delimiter>-separated string format."""
        coords = [' '.join(line.split()) for line in str(self).splitlines()]
        return delimiter.join(coords)

    @classmethod
    def from_string(cls, s, delimiter='\n', **kwargs):
        """Instantiate Molecule from <delimiter>-separated string:

            <atom1 x1 y1 z1<delimiter>atom2 x2 y2 z2<delimiter>atom3 x3 y3 z3 ...>
        """
        lines = [line for line in s.split(delimiter)]
        symbols = np.array([[atom.split()[0]] for atom in lines])
        coords = np.array([[float(c) for c in atom.split()[1:]] for atom in lines])
        n_atoms = coords.shape[0]
        comment = 'Molecule instantiated from string'

        # Instantiate the Molecule
        mol = cls(**kwargs)
        mol.coords, mol.comment, mol.symbols, mol.n_atoms = coords, comment, symbols, n_atoms
        return mol

    @staticmethod
    def load_xyzfile(xyz_file):
        """Read XYZ file and return its data."""
        with open(xyz_file) as f:
            lines = [line.strip() for line in f.readlines()]
        n_atoms = int(lines[0])
        comment = lines[1]
        symbols = np.array([[atom.split()[0]] for atom in lines[2:]])
        coords = np.array([[float(c) for c in atom.split()[1:]] for atom in lines[2:]])

        # Sanity check of XYZ file and the parsing
        assert n_atoms == coords.shape[0] == symbols.shape[0], f'Invalid XYZ file: Number of atoms does not match.'
        return n_atoms, comment, symbols, coords


class InvalidInputSection(KeyError):
    pass


class InvalidInputKeyword(KeyError):
    pass


if __name__ == '__main__':
    m = MRChemInputGenerator(hide_defaults=True)
    m.add_input_section('Molecule', 'SCF', 'WaveFunction', 'MPI')
    m.input.world_prec = 1e-4
    m.input.SCF.kain = 10
    m.input.SCF.max_iter = 100
    m.input.WaveFunction.method = 'lda'
    m.input.WaveFunction.restricted = False
    m.input.Molecule.charge = 1
    m.input.Molecule.multiplicity = 4
    m.input.MPI.numerically_exact = True

    mol = Molecule(xyz_file='NH3O.xyz')
    mol2 = Molecule.from_string(mol.to_string_format())
    assert mol == mol2, 'Molecules do not match!'
