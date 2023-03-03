import setuptools
from setuptools import find_packages
import os

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requirements = [
    # "PyAstronomy",
]

setuptools.setup(
    name="lamostabs", 
    version="0.0.1",
    author="Wang et al.",
    author_email="",
    description="Absolute flux calibration for LAMOST Low Resolution Spectral (LRS) data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wang1hm/LamostAbs",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GENERAL PUBLIC LICENSE V3",
        "Operating System :: OS Independent",
    ],
    install_requires=requirements,
    python_requires='>=3.7',
)
