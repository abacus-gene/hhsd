from setuptools import setup, find_packages
setup(
    name='hhsd',
    version='0.9.7',

    author='Daniel Kornai',
    description='Hierarchical heuristic species delimitation under the multispecies coalescent model with migration.',

    url='https://github.com/abacus-gene/hhsd',

    packages=find_packages(),

    install_requires=[
        'numpy>=1.20',
        'scipy>=1.10.1',
        'pandas>=1.5',
        'Bio>=1.5',
        ],

    python_requires = '>=3.9,<4',

    entry_points={
        'console_scripts': [
            'hhsd=hhsd.hhsd:run'
        ]
    }
)