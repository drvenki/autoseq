from setuptools import setup
from pip.req import parse_requirements

# parse_requirements() returns generator of pip.req.InstallRequirement objects
install_reqs = parse_requirements("requirements.txt", session=False)

# reqs is a list of requirement
reqs = [str(ir.req) for ir in install_reqs if ir.req is not None]

setup(name='pyautoseq',
      version='0.3.3',
      packages=['autoseq', 'autoseq.pipeline', 'autoseq.tools', 'autoseq.util'],
      install_requires=reqs,
      entry_points={
          'console_scripts': [
              'pyautoseq = autoseq.cli:main',
              'report2json = autoseq.report2json:main',
              'generate-ref = autoseq.generate_ref:main'
          ]
      }
      )
