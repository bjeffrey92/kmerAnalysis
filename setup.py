from setuptools import setup
from distutils.command.install import install as _install
import subprocess

class install(_install):
    def run(self):
        subprocess.call(['make', 'clean', '-C', 'src'])
        subprocess.call(['make', '-C', 'src'])
        _install.run(self)

setup(name='kmerAnalysis',
      version='0.1',
      description='',
      url='',
      author='Benjamin Jeffrey',
      author_email='b.jeffrey@imperial.ac.uk',
      packages=['kmerAnalysis'],
      install_requires=['pandas',
                        'numpy',
                        'bloom_filter',
                        'sklearn',
                        'matplotlib'],
      setup_requires=['wheel'],
      package_data={'kmerAnalysis': ['combinekmerdata.so']},
      cmdclass={'install': install}
)
