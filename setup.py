from setuptools import setup, find_packages
setup(
    name='hhsd',
    version='0.9.5',

    author='Daniel Kornai',
    description='Hierarchical heuristic species delimitation under the multispecies coalescent model with migration.',

    url='https://github.com/abacus-gene/hhsd',

    packages=find_packages(),

    install_requires=[
        'distinctipy',
        'ete3',
        'numpy',
        'scipy',
        'Bio',
        'pandas',
        'pyqt5',
        'lxml',
        'six',
        ],

    entry_points={
        'console_scripts': [
            'hhsd=hhsd:run'
        ]
    }
)