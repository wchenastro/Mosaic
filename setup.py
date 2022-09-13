from setuptools import setup

setup(name='mosaic',
      version='1.2.0',
      description='Code for generating point spread functions and beam tiling for radio interferometers',
      url='https://github.com/wchenastro/Mosaic',
      author='Weiwei Chen',
      author_email='wchen@mpifr-bonn.mpg.de',
      license='MIT',
      packages=['mosaic'],
      install_requires=[
          'scipy',
          'numpy',
          'matplotlib',
          'katpoint',
          'nvector == 0.7.6',
          'astropy'
      ],
      zip_safe=False)
