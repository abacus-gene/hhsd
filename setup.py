from setuptools import setup, find_packages
setup(
    name='hhsd',
    version='0.9.8',
    author='Daniel Kornai',
    description='Hierarchical heuristic species delimitation under the multispecies coalescent model with migration.',
    url='https://github.com/abacus-gene/hhsd',
    
    packages=find_packages(),

    package_data={
        'hhsd': ['bpp/*/*'],  # include all files under bpp folder and its subfolders
    },

    install_requires=[
        'numpy>=1.20',
        'scipy>=1.10.1',
        'pandas>=1.5',
        'Bio>=1.5',
    ],
    python_requires='>=3.9,<3.10',
    entry_points={
        'console_scripts': [
            'hhsd=hhsd.hhsd:run'
        ]
    }
)