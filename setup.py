# setup.py

from setuptools import setup, find_packages
from scpp.__init__ import __VERSION__


setup(
    name="scpp",
    version=__VERSION__,
    description="scpp for single cell analyse",
    author="liu chenglong",
    author_email="njlcl@outlook.com",
    install_requires=["omicverse"],
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "scpp = scpp.scpp:main",
        ],
    },
    include_package_data=True,
)
