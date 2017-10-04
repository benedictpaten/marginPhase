from version import version, required_versions
from setuptools import find_packages, setup


kwargs = dict(
    name='toil-marginPhase',
    version=version,
    description="UCSC CGL MarginPhase Toil pipeiline",
    author='UCSC Computational Genomics Lab',
    author_email='cgl-toil@googlegroups.com',
    url="https://github.com/",
    install_requires=[x + y for x, y in required_versions.iteritems()],
    tests_require=['pytest==2.8.3'],
    package_dir={'': 'src'},
    packages=find_packages('src'),
    entry_points={
        'console_scripts': ['toil-marginphase = toil_marginphase.marginphase_pipeline:main']})


setup(**kwargs)
