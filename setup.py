from setuptools import setup, find_packages

setup(
    name='mmfex',
    version='1.0.1',
    author='Roberto Cacciari, Luca De Vico',
    author_email='luca.devico@unisi.it',
    packages=find_packages(),
    url='https://github.com/m02n1ng/mmfex.git',
    description='Frenkel-Davydov exciton Hamiltonian from OpenMolcas outputs',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    python_requires='>=3.10',
    install_requires=['numpy>=2.1,<3',
                      'scipy>=1.14,<2'
                     ],
    extras_require={'test': ['pytest']}
)