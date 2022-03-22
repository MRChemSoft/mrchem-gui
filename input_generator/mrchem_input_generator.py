import copy
import json
import yaml
import requests
import numpy as np
from types import SimpleNamespace
from pprint import pformat


def couldbefloat(num):
    """Helper function"""
    try:
        float(num)
        return True
    except ValueError:
        return False


def get_type(t):
    """Helper function. Return builtin type."""
    return getattr(__builtins__, t)


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


class MRChemInputReference:
    """Docstring"""
    def __init__(self, template_file=None):
        if template_file is None:
            self.file_template = 'https://raw.githubusercontent.com/MRChemSoft/mrchem/master/python/template.yml'
        else:
            self.file_template = template_file
        self.__dict__.update(**yaml.safe_load(requests.get(self.file_template).text))
        self.defaults = self.__parse()

    def get_input_section(self, section):
        """Return dict of (key, val) pairs for passed section."""
        return {
            key: val['value'] for key, val in self.defaults[section].items()
        }

    def get_all_sections(self):
        """Return dict of (key, val) pairs for all sections."""
        top_level = {
            key: val['value'] for key, val in self.defaults.items() if 'world' in key
        }

        sections = {
            secname: {
                key: val['value'] for key, val in secval.items()
            } for secname, secval in self.defaults.items() if 'world' not in secname
        }

        return {**top_level, **sections}

    def __parse(self):
        """
        Read the MRChem input template, and generate a nested dict of all keyword
        defaults (if present), types, and docstrings
        """
        top_level = {}
        for key in self.keywords:
            name = key['name']
            docs = key['docstring']
            type = key['type']
            try:
                default = key['default']
            except KeyError:
                default = None
            top_level[name] = {}
            top_level[name]['docs'] = docs.strip()
            top_level[name]['type'] = type
            top_level[name]['value'] = default

        sections = {}
        for section in self.sections:
            secname = section['name']
            sections[secname] = {}
            for key in section['keywords']:
                keyname = key['name']
                keydocs = key['docstring']
                keytype = key['type']
                try:
                    keydefault = key['default']
                except KeyError:
                    keydefault = None
                sections[secname][keyname] = {}
                sections[secname][keyname]['docs'] = keydocs
                sections[secname][keyname]['type'] = keytype
                sections[secname][keyname]['value'] = keydefault
        return {**top_level, **sections}


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
    json_indentation = 2

    def __init__(self, hide_defaults=True, template_file=None):
        self.hide_defaults = hide_defaults
        self.reference = MRChemInputReference(template_file=template_file)
        self.input = RecursiveNamespace(**{})

        # Add the top level 'world' keywords
        for key in [key for key in self.reference.defaults.keys() if 'world' in key]:
            self.input.__setattr__(
                key,
                self.reference.defaults[key]['value']
            )

    def __str__(self):
        if not self.hide_defaults:
            return json.dumps(self.to_dict(self.input), indent=self.json_indentation)
        else:
            return json.dumps(self.__without_defaults(), indent=self.json_indentation)

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
        raise NotImplementedError('This is not implemented yet.')

    def add_input_section(self, *sections):
        """Add input section."""
        for section in sections:
            if not section in self.reference.defaults.keys():
                raise InvalidInputSection(section)
            self.input.__setattr__(
                section, RecursiveNamespace(**self.reference.get_input_section(section))
            )

    def remove_input_section(self, *sections):
        """Delete input section."""
        for section in sections:
            if not section in self.reference.defaults.keys():
                raise InvalidInputSection(section)
            self.input.__delattr__(section)

    def __without_defaults(self):
        """Return input dict without keywords with default values."""
        inp = self.to_dict(self.input)
        ref = self.reference.defaults
        new = {}

        # Only add keywords/sections with non-default values
        for section in inp:
            if isinstance(inp[section], dict):
                sub = {}
                sub[section] = {}
                for key, val in inp[section].items():
                    if val != ref[section][key]['value']:
                        sub[section][key] = val
                new.update(sub)

                if new[section] == {}:
                    del new[section]
            # Top level keywords
            else:
                if inp[section] != ref[section]['value']:
                    new[section] = inp[section]

        return new


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
    def __init__(self, charge=0, multiplicity=1, precision=10):
        self.charge = charge                                  # Total charge
        self.multiplicity = multiplicity                      # Total spin multiplicity
        self.precision = precision                            # Number of decimals in coordinate
        self.n_unpaired_electrons = self.multiplicity - 1     # Number of unpaired electrons

        self.xyz_file = None
        self.natoms = None
        self.coords = None
        self.symbols = None
        self.comment = None

    def __str__(self):
        """Pretty print coordinate section of XYZ file."""
        symbols = self.symbols.tolist()
        coords = self.coords.tolist()
        new = []
        for sym, c in zip(symbols, coords):
            x, y, z = c
            align = lambda c: f'{c:.{self.precision}f}' if c < 0 else f' {c:.{self.precision}f}'
            x = align(x)
            y = align(y)
            z = align(z)

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
    def from_xyzfile(cls, xyzfile, **kwargs):
        """Instantiate Molecule from XYZ file."""
        mol = cls(**kwargs)
        mol.xyz_file = xyzfile
        mol.n_atoms, mol.comment, mol.symbols, mol.coords = mol.load_xyzfile(xyzfile)
        return mol

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
    e = MRChemInputGenerator(hide_defaults=False)
    e.add_input_section('SCF', 'ZORA')
