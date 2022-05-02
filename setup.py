"""The setup module for Concatenator"""

# Always prefer setuptools over distutils
from setuptools import setup, find_namespace_packages
import pathlib

# Get the long description from the README file
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="concatenator",
    version="0.2.1",
    description="Performs a sequence of transformations on an input file",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/iTaxoTools/concatenator/",
    author="Vladimir Kharchev",
    package_dir={"": "src"},
    packages=find_namespace_packages(
        include=("itaxotools*",),
        where="src",
    ),
    python_requires=">=3.8.6, <4",
    install_requires=[
        "zipp>=3.6.0",  # BUGFIX: backport from Python 3.9.1
        "pandas>=1.3.0",
        "regex>=2021.8.28",
        "networkx>=2.8",
        "itaxotools-common==0.2.2",
        "DNAconvert>=0.2.0",
    ],
    extras_require={
        "dev": [
            "pyinstaller>=4.5.1",
            "pytest>=6.2.5",
        ],
    },
    entry_points={
        "console_scripts": [
            "concatenator=itaxotools.concatenator:main",
        ]
    },
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3 :: Only",
    ],
    include_package_data=True,
)
