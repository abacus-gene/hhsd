import os
import stat
from setuptools import setup, find_packages
from setuptools.command.install import install

class PostInstallCommand(install):
    """Custom install command to handle post-install executable setting for bpp"""
    def run(self):
        install.run(self)

        base_dir = os.path.join(self.install_lib, 'hhsd', 'bpp')
        for root, dirs, files in os.walk(base_dir):
            for file in files:
                file_path = os.path.join(root, file)
                try:
                    st = os.stat(file_path)
                    os.chmod(file_path, st.st_mode | stat.S_IEXEC)
                    print(f'Set executable permission for {file_path}')
                except Exception as e:
                    print(f'Failed to set executable permission for {file_path}: {e}')


setup(
    name='hhsd',
    version='1.1.0',
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
    python_requires='>=3.9,<=3.14',
    entry_points={
        'console_scripts': [
            'hhsd=hhsd.hhsd:run'
        ]
    },
    cmdclass={'install': PostInstallCommand},
)
