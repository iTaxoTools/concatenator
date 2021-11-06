"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_namespace_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='concatenator',

    version='0.2.dev0',

    description='concatenator description',

    long_description=long_description,

    long_description_content_type='text/markdown',

    url='https://github.com/iTaxoTools/concatenator/',

    author='Vladimir Kharchev',

    # Classifiers help users find your project by categorizing it.
    #
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[  # Optional
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],

    package_dir={'': 'src'},

    packages=find_namespace_packages(
        include=('itaxotools*',),
        where='src',
    ),

    python_requires='>=3.8, <4',

    install_requires=[
        'zipp>=3.6.0',  # BUGFIX: backport from Python 3.9.1
        'pandas>=1.3.0',
        'regex>=2021.8.28',
        'itaxotools-common==0.2.dev2',
        'DNAconvert>=0.1.dev1',
    ],

    extras_require={
        'dev': [
            'pyinstaller>=4.5.1',
            'pytest>=6.2.5',
        ],
    },

    # Include all data from MANIFEST.in
    include_package_data=True,

    entry_points={
        'console_scripts': [
            'concatenator=itaxotools.concatenator:main',
        ],
        'pyinstaller40': [
            'hook-dirs = itaxotools.__pyinstaller:get_hook_dirs',
            'tests = itaxotools.__pyinstaller:get_pyinstaller_tests'
        ]
    },

)
