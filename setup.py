from setuptools import setup, find_packages
setup(
    name='hhsd',
    version='0.9.5',

    author='Daniel Kornai',
    description='Hierarchical heuristic species delimitation under the multispecies coalescent model with migration.',

    url='https://github.com/abacus-gene/hhsd',

    packages=find_packages(),

    install_requires=[
        'distinctipy>=1.2.1',
        'ete3>=3.1.2',
        'numpy>=1.20',
        'scipy>=1.10.1',
        'Bio>=1.5',
        'pandas>=1.5',
        'pyqt5>=5.15.9',
        'lxml>=4.9.2',
        'six>=1.16',
        ],

    python_requires = '>=3.9,<4',

    entry_points={
        'console_scripts': [
            'hhsd=hhsd:run'
        ]
    }
)
