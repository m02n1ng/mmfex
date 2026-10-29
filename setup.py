from setuptools import setup, find_packages

setup(
    name='mmfex',
    version=1.0,
    author='Roberto Cacciari',
    author_email='cacciari.roberto@phd.unipd.it',
    packages=find_packages(),
    url='https://github.com/m02n1ng/mmfex.git',
    install_requires=['numpy==2.1.3',
                      'scipy==1.14.1'
                     ],
    extras_require={'test': ['pytest']}
)