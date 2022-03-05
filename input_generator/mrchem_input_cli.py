#!/usr/bin/env python

import argparse
from mrchem_input_generator import MRChemInputGenerator
from mrchem_input_templates import EnergyCalculation, ElectricResponseCalculation, MagneticResponseCalculation

CURRENT_INPUT_TYPES = {
    'energy': {
        'class': EnergyCalculation,
        'default_filename': 'energy.inp'
    },
    'electric': {
        'class': ElectricResponseCalculation,
        'default_filename': 'electric.inp'
    },
    'magnetic': {
        'class': MagneticResponseCalculation,
        'default_filename': 'magnetic.inp'
    }
}


def cli():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     prog='mrchem_input_templates.py',
                                     description='',
                                     epilog='')
    subparsers = parser.add_subparsers(dest='command')

    new = subparsers.add_parser('new')
    new.add_argument('input_type',
                     type=str,
                     default='energy',
                     choices=CURRENT_INPUT_TYPES.keys(),
                     help='What type of job do you want?')
    new.add_argument('-f', '--filename',
                     type=str,
                     help='Filename for generated input file.')
    new.add_argument('-d', '--with_defaults',
                     action='store_true',
                     help='Include keywords with default values.')
    new.add_argument('-xyz', '--xyzfile', help='Path to XYZ file.')
    new.add_argument('-p', '--print', action='store_true', help='Write to stdout instead of to file.')

    add = subparsers.add_parser('add')
    add.add_argument('input_file',
                     help='Path to input file to be modified.')
    add.add_argument('-a', '--add', nargs=3, help='3-tuple indicating which keyword to add/update <sec key val>')
    add.add_argument('-xyz', '--xyzfile', help='Path to XYZ file.')

    return parser


if __name__ == '__main__':
    parser = cli()
    args = parser.parse_args()

    # Make new input file
    if args.command == 'new':
        if args.filename is None:
            args.fname = CURRENT_INPUT_TYPES[args.input_type]['default_filename']

        mrc = CURRENT_INPUT_TYPES[args.input_type]['class'](hide_defaults=not args.with_defaults,
                                                            xyzfile=args.xyzfile)
        if args.print:
            print(mrc)
        else:
            with open(args.fname, 'w') as f:
                f.write(str(mrc))

    # Modify existing input file
    elif args.command == 'add':
        mrc = MRChemInputGenerator.from_json(args.input_file)
        if args.add is not None:
            section, key, value = args.add
            mrc.set_keyword(section, key, value)

        # Add coordinates
        if args.xyzfile is not None:
            coords, _, _ = mrc._load_xyzfile(args.xyzfile)
            mrc.set_keyword('Molecule', 'coords', coords)

        # Write updated input to file
        with open(args.input_file, 'w') as f:
            f.write(str(mrc))