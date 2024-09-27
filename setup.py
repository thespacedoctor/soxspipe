from setuptools import setup, find_packages
import os

moduleDirectory = os.path.dirname(os.path.realpath(__file__))
exec(open(moduleDirectory + "/soxspipe/__version__.py").read())


def readme():
    with open(moduleDirectory + '/README.md') as f:
        return f.read()


install_requires = [
    'pyyaml==5.3.1',
    'soxspipe',
    'fundamentals',
    'astropy',
    'ccdproc',
    'docopt',
    'photutils',
    'matplotlib',
    'numpy',
    'unicodecsv',
    'pandas',
    'tabulate',
    'bottleneck',
    'multiprocess',
    'jinja2',
    'specutils'
]

# READ THE DOCS SERVERS
exists = os.path.exists("/home/docs/")
if exists:
    install_requires = ['fundamentals']


setup(name="soxspipe",
      version=__version__,
      description="A python package and command-line tools to The data-reduction pipeline for the SOXS instrument",
      long_description=readme(),
      long_description_content_type='text/markdown',
      classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 3.9',
          'Topic :: Utilities',
      ],
      keywords=['soxs, eso, data, pipeline, spectra'],
      url='https://github.com/thespacedoctor/soxspipe',
      download_url='https://github.com/thespacedoctor/soxspipe/archive/v%(__version__)s.zip' % locals(
      ),
      author='David Young',
      author_email='davidrobertyoung@gmail.com',
      license='GPLv3',
      packages=find_packages(exclude=["*tests*"]),
      include_package_data=True,
      install_requires=install_requires,
      test_suite='nose2.collector.collector',
      tests_require=['nose2', 'cov-core'],
      entry_points={
          'console_scripts': ['soxspipe=soxspipe.cl_utils:main'],
      },
      zip_safe=False)
