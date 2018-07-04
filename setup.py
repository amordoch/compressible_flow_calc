from setuptools import setup
from compressible_flow_calc import __CFCVERSION__

setup(
    name='compressible_flow_calc',
    version=__CFCVERSION__,
    packages=['compressible_flow_calc'],
    install_requires=['numpy'],
    url='',
    license='GPLv3',
    author='Ariel Mordoch',
    author_email='arielmordoch@gmail.com',
    description='Compressible flow calculator for aerodynamics'
)
