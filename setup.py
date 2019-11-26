from setuptools import setup

setup(name='kmerAnalysis',
      version='0.1',
      description='The funniest joke in the world',
      url='',
      author='Benjamin Jeffrey',
      author_email='b.jeffrey@imperial.ac.uk',
      license='MIT',
      packages=['kmerAnalysis'],
      install_requires=['pandas',
                        'numpy',
                        'bloom_filter',
                        'sklearn',
                        'matplotlib']
)
