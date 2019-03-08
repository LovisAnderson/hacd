from setuptools import setup

config = {
    'name': 'hacd',
    'description': 'Package to calculate a hierarchical approximative convex decomposition'
                   ' of a union of polytopes.',
    'author': 'Lovis Anderson',
    'author_email': 'lovisanderson@posteo.net',
    'license': 'GPL',
    'version': '0.1',
    'packages': ['hacd'],
    'package_dir': {'hacd': 'hacd'},
    'install_requires': [
        'enum34',
        'sweepvolume',
        'argparse',
    ],
}

setup(**config)