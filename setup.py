from setuptools import setup

exec(open('main/version.py').read())

setup(
    name= 'DeNovoFilter',
    version=__version__,
    description='A tool for identifying potential de novo variants',
    url='https://github.com/RahmanTeamDevelopment/DeNovoFilter',
    author='RahmanTeam',
    author_email='rahmanlab@icr.ac.uk',
    license='MIT',
    packages=['main'],
    scripts=['bin/DeNovoFilter.py'],
    zip_safe=False
)
