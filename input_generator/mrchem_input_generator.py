import copy
import json
import yaml
import requests
from types import SimpleNamespace


class PrettyDict(dict):
    """Subclass of dict to enforce default pretty printing.
    Perhaps best to avoid this and let the user pretty print if needed."""
    def __init__(self, **d):
        super().__init__(**d)

    def __str__(self):
        return json.dumps(self, indent=4)


class RecursiveNamespace(SimpleNamespace):
    @staticmethod
    def map_to_ns(entry):
        """Map [dict] to [SimpleNamespace]"""
        if isinstance(entry, dict):
            return RecursiveNamespace(**entry)

    @staticmethod
    def map_to_dict(entry):
        """Map [SimpleNamespaces] to [dict]"""
        if isinstance(entry, SimpleNamespace):
            return entry.__dict__

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        for key, val in kwargs.items():
            if isinstance(val, dict):
                setattr(self, key, RecursiveNamespace(**val))
            elif isinstance(val, list):
                setattr(self, key, val)


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
            null,
            null,
            null
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
    _remote_template = 'https://raw.githubusercontent.com/MRChemSoft/mrchem/release/1.0/python/template.yml'
    _json_indentation = 2

    def __init__(self, hide_defaults=True):
        self.hide_defaults = hide_defaults
        self.template = yaml.safe_load(requests.get(self._remote_template).text)

        self.input = RecursiveNamespace(**{})
        self._defaults = RecursiveNamespace(**self._parse_tmplt())

        # Add the top level 'world' keywords
        for key in [key for key in self.to_dict(self._defaults).keys() if 'world' in key]:
            self.input.__setattr__(key, self._defaults.__getattribute__(key))

    def __str__(self):
        """Pretty string format for printing and/or file writing."""
        return json.dumps(self.generate_without_defaults(), indent=self._json_indentation)

    def _parse_tmplt(self):
        """
        Read the MRChem input template, and generate a nested dict of all keywords
        and their defaults (if present), which can easily be turned in to a recursive namespace
        (i.e. to get dotted access to all keys).
        :return: dict
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
        """Add input section to object and list of attributes."""
        for section in sections:
            self.input.__setattr__(
                section, RecursiveNamespace(**self._defaults.__getattribute__(section).__dict__)
            )

    def remove_input_section(self, *sections):
        """Delete input section from object and list of attributes."""
        for section in sections:
            self.input.__delattr__(section)

    @staticmethod
    def to_dict(ns):
        """Recursively convert SimpleNamespaces to dictionaries."""
        d = {}
        for key, val in ns.__dict__.items():
            if isinstance(val, SimpleNamespace):
                d[key] = val.__dict__
            else:
                d[key] = val
        return PrettyDict(**d)

    def generate_with_defaults(self):
        """Return pretty input with default keywords."""
        return self.to_dict(self.input)

    def generate_without_defaults(self):
        """Return pretty input without default keywords."""
        i = self.to_dict(self.input)
        d = self.to_dict(self._defaults)

        keys_to_delete = []
        filtered = copy.deepcopy(i)
        for section in i.keys():
            if isinstance(i[section], dict):
                for key, keyval in i[section].items():
                    val = i[section][key]
                    ref = d[section][key]
                    if val == ref:
                        keys_to_delete.append((section, key))

        for section, key in keys_to_delete:
            del filtered[section][key]
            if filtered[section] == {}:
                del filtered[section]

        return filtered


if __name__ == '__main__':
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
