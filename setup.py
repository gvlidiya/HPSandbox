
#!/usr/bin/env python
from setuptools import setup, find_packages
from os import path
from io import open

here = path.abspath(path.dirname(__file__))
requirements_path = path.join(here, "requirements.txt")

# Read requirements.txt without relying on pip internals (works in build isolation)
requirements_path = path.join(here, "requirements.txt")
with open(requirements_path, "r", encoding="utf-8") as f:
    requirements = [
        line.strip()
        for line in f
        if line.strip() and not line.strip().startswith("#")
    ]

# Get the long description from the README file
with open(path.join(here, "README.md"), "r", encoding="utf-8") as f:
    long_description = f.read()


setup(
        name='hpsandbox',
        version='1.0.0b',
        #scripts=['hpsandbox'] ,
        author="Vincent Voelz",
        author_email="vvoelz@gmail.com",
        description="Python package to explore the 2D HP model of Chan and Dill.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/vvoelz/HPSandbox",
        packages=find_packages(),
        install_requires=requirements,
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ]
)





