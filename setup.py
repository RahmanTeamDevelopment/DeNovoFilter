from setuptools import setup

setup(
    name= 'DeNovoFilter',
    version = '0.1.0',
    description = 'A tool for identifying potential de novo variants',
    url = 'https://github.com/RahmanTeamDevelopment/DeNovoFilter',
    author = 'RahmanTeam',
    author_email = 'rahmanlab@icr.ac.uk',
    license = 'MIT',
    packages=['denovo_'],
    scripts=['bin/DeNovoFilter.py','bin/denovo'],
    zip_safe=False
)
