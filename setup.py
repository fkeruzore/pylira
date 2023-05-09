from setuptools import setup

setup(
    name="pylira",
    version="0.1.2",
    description="pylira: Python wrapper for the LIRA R library",
    author="Florian Keruzore",
    author_email="fkeruzore@anl.gov",
    packages=["pylira"],
    install_requires=[
        "numpy",
        "matplotlib",
        "astropy",
        "pandas",
        "scipy",
        "chainconsumer",
        "rpy2",
    ],
    include_package_data=True,
    package_data={"": ["./r_script.R"]},
)
